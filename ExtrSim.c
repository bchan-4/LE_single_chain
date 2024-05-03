#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#define LAMMPS_LIB_MPI
#include "library.h" /* this is a LAMMPS include file */

double UBond(double K, double R0, double currR)
{ // function to calculate bond energy
    double U;
    U = -0.5 * K * R0 * R0 * log(1 - (currR / R0) * (currR / R0)) + 4 * ((1 / pow(currR, 12)) - (1 / pow(currR, 6))) + 1;
    return U;
}

double getdist(double *coord1, double *coord2)
{ // Get distance between two beads, without periodic boundaries
    double distvect[3];
    for (int xx = 0; xx < 3; xx++)
    {
        distvect[xx] = coord1[xx] - coord2[xx];
    }
    double dist = sqrt(distvect[0] * distvect[0] + distvect[1] * distvect[1] + distvect[2] * distvect[2]);
    return dist;
}

int CurrNumCoh(int *cohmat, int maxnum)
{ // function to get number of bound cohesins. Goes through FIRST HALF of the matrix and sees where >0
    int numbCoh = 0;
    for (int kk = 0; kk < maxnum; kk++)
    {
        if (cohmat[kk] > 0)
        {
            numbCoh = numbCoh + 1;
        }
    }
    return numbCoh;
}

int *getwherebound(int *indbound, int *cohlifes, int maxnum)
{ // function to get indices where cohesins are bound. Can also be used to see which indices of newbonds/delbonds array are filled.
    int boundindex = 0;
    for (int kk = 0; kk < maxnum; kk++)
    {
        if (cohlifes[kk] > 0)
        {
            indbound[boundindex] = kk + 1; // Need to add 1 because LAMMPS starts indexing from 1
            boundindex = boundindex + 1;
        }
    }
    return indbound;
}

double *getcoords(double *thiscoord, int beadindex, double *allcoords)
{ // get coordinates of specific bead. Bead index is in LAMMPS indexing (1 to N)
    for (int ii = 0; ii < 3; ii++)
    {
        thiscoord[ii] = allcoords[(beadindex - 1) * 3 + ii];
    }
    return thiscoord;
}

double getE(double *AllCoords, int *CohInd, int NumBeads, double K, double R0)
{ // gets energy of a cohesin bond.
    double thisE = 0;
    double *coordsA;
    double *coordsB;
    coordsA = (double *)malloc(3 * sizeof(double));
    coordsB = (double *)malloc(3 * sizeof(double));
    // get coordinates of trial beads
    coordsA = getcoords(coordsA, CohInd[0], AllCoords);
    coordsB = getcoords(coordsB, CohInd[1], AllCoords);
    double Sep1 = getdist(coordsA, coordsB);
    if (Sep1 > (R0 - 0.005))
    {                    // if distance too large
        thisE = -100000; // Reject
    }
    else if (Sep1 >= 0)
    {
        thisE = UBond(K, R0, Sep1);
    }
    else
    {
        thisE = -100000;
    }
    free(coordsA);
    free(coordsB);
    return thisE;
}

int *getNewBond(int *NewLoc, double JumpProb, int *OldBond, int *CohLocs, double *AllCoords, int *yesstop, int NumBeads, int MaxCoh, double K, double R0, int yesring)
{ // function to decide indices of new cohesin bond
    double ProbFactor;
    double ThisRand = (double)rand() / (RAND_MAX); // to decide which domains to move
    int movevect[2];                               // movevect[i]=-1,1 for moving to end or to center
    int oppyesstop[2];                             // not yes stop (ie 0 becomes 1, and 1 becomes 0)
    for (int rr = 0; rr < 2; rr++)
    {
        if (yesstop[rr] == 0)
        {
            oppyesstop[rr] = 1;
        }
        else
        {
            oppyesstop[rr] = 0;
        }
    }
    // always active for now
    movevect[0] = -1;
    movevect[1] = 1;

    // which side to try and move
    int whichsides[2];
    if (ThisRand <= 0.25)
    {
        whichsides[0] = 0;
        whichsides[1] = 1;
    }
    else if (ThisRand <= 0.5)
    {
        whichsides[0] = 1;
        whichsides[1] = 0;
    }
    else if (ThisRand <= 0.75)
    {
        whichsides[0] = 1;
        whichsides[1] = 1;
    }
    else
    {
        whichsides[0] = 0;
        whichsides[1] = 0;
    }

    // Moving steps
    if (yesstop[0] == 1 && yesstop[1] == 1)
    { // if currently stopped;
        NewLoc[0] = OldBond[0];
        NewLoc[1] = OldBond[1];
    }
    else
    {                     // otherwise, try to move the bond
        int Move1temp[2]; // initial trial for moving bond
        for (int ii = 0; ii < 2; ii++)
        {
            Move1temp[ii] = OldBond[ii] + whichsides[ii] * movevect[ii] * oppyesstop[ii];
        }
        if (yesring == 0)
        {
            if (Move1temp[0] < 1 || Move1temp[1] > NumBeads || Move1temp[0] > NumBeads || Move1temp[1] < 1)
            {
                // If i try to move past ends, just fall off
                NewLoc[0] = 0;
                NewLoc[1] = 0;
                ProbFactor = 1;
                return NewLoc;
            }
        }
        else // if ring
        {    // wrap around
            if (Move1temp[0] == 0)
            {
                Move1temp[0] = NumBeads;
            }
            if (Move1temp[1] > NumBeads)
            {
                Move1temp[1] = 1;
            }

            // check if went all the way around
            if (Move1temp[0] == OldBond[1] || Move1temp[1] == OldBond[0])
            {
                NewLoc[0] = 0;
                NewLoc[1] = 0;
                ProbFactor = 1;
                return NewLoc;
            }
        }
        // Check if need to stall
        int yesmove[2];
        yesmove[0] = 1;
        yesmove[1] = 1;
        for (int ww = 0; ww < 2; ww++)
        { // for each domain
            if (Move1temp[ww] * whichsides[ww] * oppyesstop[ww] > 0)
            {                       // if try to move
                int isoccupied = 0; // flag for if trial bead is occupied
                for (int vv = 0; vv < MaxCoh; vv++)
                { // go through CohLocs and see if trial bead is already occupied
                    if (CohLocs[vv] == Move1temp[ww] || CohLocs[vv + MaxCoh] == Move1temp[ww])
                    {
                        // if occupied, set isoccupied=1 and break
                        isoccupied = 1;
                        break;
                    }
                }
                if (isoccupied > 0)
                { // if something already there
                    double jumprand = (double)rand() / (RAND_MAX);
                    if (jumprand <= JumpProb)
                    { // only move based on jump probability
                        yesmove[ww] = 1;
                    }
                    else
                    {
                        yesmove[ww] = 0;
                    }
                }
            }
        } // finish stalling calc

        // Adjust Move location with stalling prob
        int Move1[2];
        for (int ff = 0; ff < 2; ff++)
        {
            if (yesmove[ff] == 0)
            {
                Move1[ff] = OldBond[ff];
            }
            else
            {
                Move1[ff] = Move1temp[ff];
            }
        }
        // if both elements of Move1 are same bead, reject and keep old bond. Need for passive.
        if (Move1[0] == Move1[1])
        {
            NewLoc[0] = OldBond[0];
            NewLoc[1] = OldBond[1];
            return NewLoc; // exit function early
        }

        // get probability of moving
        double currE = getE(AllCoords, OldBond, NumBeads, K, R0);
        double moveE = getE(AllCoords, Move1, NumBeads, K, R0);

        // active
        if (moveE < 0)
        {
            ProbFactor = 0;
        }
        else if (moveE >= 0)
        { // move it UNLESS I ran into an issue like bond too long
            ProbFactor = 1;
        }
        else
        {
            ProbFactor = 0;
        }

        // Actually set NewLocs
        double moverand = (double)rand() / (RAND_MAX);
        if (moverand <= ProbFactor)
        {
            NewLoc[0] = Move1[0];
            NewLoc[1] = Move1[1];
        }
        else
        {
            NewLoc[0] = OldBond[0];
            NewLoc[1] = OldBond[1];
        }
    }
    return NewLoc;
}

int *getcontacts(int *matrixforcont, double *allcoords, int Nb, double contactrad)
{
    // Add to contact matrix (matrixforcont) given coordinates (allcoords), number of beads (Nb), within contact radius (contactrad)
    // matrixforcont is in row major order
    double *iicoords;
    double *jjcoords;
    iicoords = (double *)malloc(3 * sizeof(double));
    jjcoords = (double *)malloc(3 * sizeof(double));
    double thisdist;
    int thisentry;
    for (int ii = 1; ii <= Nb; ii++)
    { // for each bead. Bead index starts from 1
        iicoords = getcoords(iicoords, ii, allcoords);
        for (int jj = ii; jj <= Nb; jj++)
        { // for other beads from ii through end
            jjcoords = getcoords(jjcoords, jj, allcoords);
            thisdist = getdist(iicoords, jjcoords);
            thisentry = matrixforcont[(ii - 1) * Nb + (jj - 1)];
            if (thisdist <= contactrad)
            {
                matrixforcont[(ii - 1) * Nb + (jj - 1)] = thisentry + 1;
            }
        }
    }
    free(iicoords);
    free(jjcoords);
    return matrixforcont;
}

double *getR2mat(double *matrixforR2, double *allcoords, int Nb)
{
    // Add to internal distances r2 matrix (matrixforR2) given coordinates (allcoords), number of beads (Nb)
    // matrixforR2 is in row major order
    double *iicoords;
    double *jjcoords;
    iicoords = (double *)malloc(3 * sizeof(double));
    jjcoords = (double *)malloc(3 * sizeof(double));
    double thisdist;
    int tempind;
    for (int ii = 1; ii <= Nb; ii++)
    { // for each bead. Bead index starts from 1
        iicoords = getcoords(iicoords, ii, allcoords);
        for (int jj = ii; jj <= Nb; jj++)
        { // for other beads from ii through end
            tempind = (ii - 1) * Nb + (jj - 1);
            jjcoords = getcoords(jjcoords, jj, allcoords);
            thisdist = getdist(iicoords, jjcoords);
            matrixforR2[tempind] = matrixforR2[tempind] + thisdist * thisdist;
        }
    }
    free(iicoords);
    free(jjcoords);
    return matrixforR2;
}

int *updatestopflags(int *stopflags, int statleft, int statright, int *newLocs, int numstoppairs, int *leftstoplocs, int *rightstoplocs, int cohind, int maxCoh)
{ // update anchor flags (1 is occupied, 0 is unoccupied)
    for (int hh = 0; hh < numstoppairs; hh++)
    { // for each pair of stop signs, left
        // active
        if (leftstoplocs[hh] == newLocs[0])
        { // check if left side stopped
            stopflags[cohind - 1] = 1;
        }
        if (rightstoplocs[hh] == newLocs[1])
        { // check if right side stopped
            stopflags[cohind - 1 + maxCoh] = 1;
        }
    }
    return stopflags;
}

int *generatepermute(int *finallist, int N)
{ // generate array with random permutation from 1:N
    for (int w = 0; w < N; w++)
    {
        finallist[w] = w + 1;
    }
    for (int w = 0; w < N; w++)
    {
        int j, k;
        j = rand() % (N - w) + w; // random location
        k = finallist[j];
        finallist[j] = finallist[w];
        finallist[w] = k; // Swap w and j
    }
    return finallist;
}

int main(int narg, char **argv)
{
    // Setup MPI
    int myrank, nprocs, provided;
    MPI_Init_thread(&narg, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    fprintf(stdout, "Hello from the beginning, rank %d\n", myrank);
    // decide which rank will handle status updates, contacts and internal distance calculations
    int statusrank = 0;
    int contactrank = 0;
    int R2rank = 0;
    if (nprocs > 2)
    {
        contactrank = 1;
        R2rank = 2;
    }
    else if (nprocs > 1)
    {
        contactrank = 1;
        R2rank = 1;
    }

    // Parse through inputs to get simulation parameters
    int NumBeads, NumBlocks, initTau, loadsite, contactr2freq, yesring, maxCoh;
    int numstoppairs = 0;
    int randseed;
    double jumpprob, K, R0, tint, TauPerDump, rc, MDperCheck, onprob, offprob;
    char *WorkPath;
    char *outputfilenames;
    char *initdataFileName;
    char *allstopName;
    FILE *allstopFile;

    maxCoh = atoi(argv[1]);         // maximum number of bound extruders allowed
    jumpprob = atof(argv[2]);       // passing probability
    NumBeads = atoi(argv[3]);       // Number of beads in system
    NumBlocks = atoi(argv[4]);      // Number of MC steps
    outputfilenames = argv[5];      // prefix for output files
    K = atof(argv[6]);              // spring constant for cohesin bond
    R0 = atof(argv[7]);             // R0 for cohesin bond
    MDperCheck = atof(argv[8]);     // number of MD tau per MC step
    TauPerDump = atof(argv[9]);     // Number of tau per dumpstep for beadcoordinates.
    WorkPath = argv[10];            // where output files should be stored
    initdataFileName = argv[11];    // name of initial data file
    initTau = atoi(argv[12]);       // Init MD tau before MC
    tint = atof(argv[13]);          // integration timestep (usually 0.01)
    loadsite = atoi(argv[14]);      // 0 if nonspecific loading, otherwise the site of loading.
    rc = atof(argv[15]);            // contact distance cutoff
    contactr2freq = atoi(argv[16]); // every contactr2freq steps, print contact matrix and cumulative R2 matrix
    yesring = atoi(argv[17]);
    randseed = atoi(argv[18]);
    onprob = atof(argv[19]);
    offprob = atof(argv[20]);
    int initdt = (int)((double)initTau / tint); //

    // print inputs just to make sure
    char statusfilename[500];
    FILE *statusFile;
    sprintf(statusfilename, "%sStatus_%s.txt", WorkPath, outputfilenames);
    if (myrank == statusrank)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 423\n");
            return 1;
        }
        fprintf(statusFile, "jumpprob=%f\n", jumpprob);
        fprintf(statusFile, "onprob=%f\n", onprob);
        fprintf(statusFile, "offprob=%f\n", offprob);
        fprintf(statusFile, "NumBeads=%d\n", NumBeads);
        fprintf(statusFile, "NumBlocks=%d\n", NumBlocks);
        fprintf(statusFile, "MDperCheck=%.3f, TauPerDump=%.3f\n", MDperCheck, TauPerDump);
        fprintf(statusFile, "K=%f, R0=%f\n", K, R0);
        fprintf(statusFile, "initTau=%d, tint=%f\n", initTau, tint);
        fprintf(statusFile, "rc=%f, contactr2freq=%d\n", rc, contactr2freq);
        fprintf(statusFile, "yesring=%d\n", yesring);
        fclose(statusFile);
    }

    // Set stop sign locations
    int *leftstoplocs;
    int *rightstoplocs;
    int withstop;
    withstop = 0;
    if (narg > 21)
    { // If stop signs provided
        withstop = 1;
        allstopName = argv[21];        // name of file with stop signs.
        numstoppairs = atoi(argv[22]); // number of stop sign pairs
        // read stop sign locations, store left and right locs in corresponding arrays
        leftstoplocs = (int *)malloc(sizeof(int) * numstoppairs);
        rightstoplocs = (int *)malloc(sizeof(int) * numstoppairs);
        allstopFile = fopen(allstopName, "r");
        if (allstopFile == NULL)
        {
            fprintf(stdout, "allstop failed line 455\n");
            return 1;
        }
        int temploc1, temploc2;
        for (int i = 0; i < numstoppairs; i++)
        { // for all stop pairs, read locs
            fscanf(allstopFile, "%d %d\n", &temploc1, &temploc2);
            leftstoplocs[i] = temploc1;
            rightstoplocs[i] = temploc2;
        }
        fclose(allstopFile);
        // make sure I read stop signs right
        if (myrank == statusrank)
        {
            statusFile = fopen(statusfilename, "a");
            if (statusFile == NULL)
            {
                fprintf(stdout, "statusfile failed line 472\n");
                return 1;
            }
            fprintf(statusFile, "First stop sign pair at %d and %d\n", leftstoplocs[0], rightstoplocs[0]);
            fprintf(statusFile, "Num stop pairs = %d\n", numstoppairs);
            fclose(statusFile);
        }
    }
    else
    { // if no stop signs, set to -1
        leftstoplocs = (int *)malloc(sizeof(int));
        rightstoplocs = (int *)malloc(sizeof(int));
        leftstoplocs[0] = -1;
        rightstoplocs[0] = -1;
    }

    // Finished parsing inputs, now do remaining setup of system
    srand((unsigned)arraynum + 1); // random number generator seed

    int *cohlocs; // holds cohesin locations. row major order
    cohlocs = (int *)malloc(sizeof(int) * maxCoh * 2);
    // initialize cohlocs
    for (int ii = 0; ii < 2 * maxCoh; ii++)
    {
        cohlocs[ii] = 0;
    }
    int *CohLifetime;
    CohLifetime = (int *)malloc(sizeof(int) * maxCoh); // holds coh lifetimes.
    for (int ii = 0; ii < maxCoh; ii++)
    {
        CohLifetime[ii] = 0;
    }
    int numCoh = 0; // number of cohesins curently bound
    int *stopflags;
    stopflags = (int *)malloc(sizeof(int) * maxCoh * 2); // 1 if stopped, 0 if not stopped. row major order
    // initialize
    for (int jj = 0; jj < 2 * maxCoh; jj++)
    {
        stopflags[jj] = 0;
    }

    int dtperdump = (int)(TauPerDump / tint); // number of integration steps per dump
    int dtperMC = (int)(MDperCheck / tint);   // number of integration steps per MC step
    // character arrays to hold LAMMPS commands
    char runcommand[100];
    char dumpcommand[500];
    char dumpmodcommand[100];
    char undumpcommand[100];
    // print out time information
    if (myrank == statusrank)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 526\n");
            return 1;
        }
        fprintf(statusFile, "dtperdump=%d, dtperMC=%d\n", dtperdump, dtperMC);
        fclose(statusFile);
    }

    // command to run LAMMPS
    sprintf(runcommand, "run %d", dtperMC);
    // command to dump coordinates
    sprintf(dumpcommand, "dump 2 all custom %d %scoordsLEUnwrap_%s.txt id type xu yu zu", dtperdump, WorkPath, outputfilenames);
    sprintf(dumpmodcommand, "dump_modify 2 sort id append yes");

    // Set up output files
    FILE *FinalLifeFile; // Final residence times of cohesins (before unbinding)
    FILE *NumCohFile;    // number of cohesins bound on the chain
    FILE *AccFile;       // Tracks acceptance of cohesin trial moves
    FILE *BindFile;      // MC step and indices when/where cohesin bound.
    char finallifename[500];
    char numcohname[500];
    char accfilename[500];
    char bindfilename[500];
    if (myrank == 0)
    {
        sprintf(finallifename, "%sFinalLife_%s.txt", WorkPath, outputfilenames);
        sprintf(numcohname, "%sNumCoh_%s.txt", WorkPath, outputfilenames);
        sprintf(accfilename, "%sAcc_%s.txt", WorkPath, outputfilenames);
        sprintf(bindfilename, "%sBind_%s.txt", WorkPath, outputfilenames);
    }
    // file to hold contacts per chain. intermittently updated
    FILE *ContactFile;
    char contactname[500];
    if (myrank == contactrank)
    {
        sprintf(contactname, "%sContacts_%s.txt", WorkPath, outputfilenames);
    }
    // file to hold internal distances per chain. intermittently updated
    FILE *R2File;
    char R2name[500];
    if (myrank == R2rank)
    {
        sprintf(R2name, "%sR2_%s.txt", WorkPath, outputfilenames);
    }
    int ffcc = 0;
    MPI_Barrier(MPI_COMM_WORLD); // Hold everything at end of setup
    if (myrank == statusrank)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 576\n");
            return 1;
        }
        fprintf(statusFile, "Finished initial setup\n");
        ffcc = fclose(statusFile);
        if (ffcc != 0)
        {
            fprintf(stdout, "fclose failed 583");
            return (1);
        }
    }

    // Start LAMMPS section
    // Command line arguments
    const char *lmpargv[] = {"liblammps", "-screen", "none", "-log", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(const char *);

    // set up LAMMPS instance
    void *lmp = NULL;
    lmp = lammps_open(lmpargc, (char **)lmpargv, MPI_COMM_WORLD, NULL);
    if (lmp == NULL)
    {
        if (myrank == 0)
        {
            statusFile = fopen(statusfilename, "a");
            if (statusFile == NULL)
            {
                fprintf(stdout, "statusfile failed line 603\n");
                return 1;
            }
            fprintf(statusFile, "LAMMPS initialization failed");
            ffcc = fclose(statusFile);
            if (ffcc != 0)
            {
                fprintf(stdout, "fclose failed 610\n");
                return (1);
            }
        }
        lammps_close(lmp);
        return 1;
    }
    else
    {
        if (myrank == 0)
        {
            statusFile = fopen(statusfilename, "a");
            if (statusFile == NULL)
            {
                fprintf(stdout, "statusfile failed line 624\n");
                return 1;
            }
            fprintf(statusFile, "LAMMPS initialized\n");
            ffcc = fclose(statusFile);
            if (ffcc != 0)
            {
                fprintf(stdout, "fclose failed 631");
                return (1);
            }
        }
    }
    // set up variable commands
    char vFileNames[100]; // Suffix for output files
    char vWorkPath[100];  // path of directory where outputs are printed
    char vK[100];         // spring constant of cohesin bind
    char vR0[100];        // R0 for cohesin bind
    char vInitFile[100];  // initial data file to set bead coordinates
    char vinitsteps[100]; // number of initial integration time steps before extrusion
    char vstepsize[100];  // timestep (dt)
    char vxdiff[100];     // used for end to end distance
    char vydiff[100];
    char vzdiff[100];
    char vRg[100];       // Radius of gyration
    char vLangSeed[200]; // random number for langevin dynamics

    int LangSeedRand = 0;
    while (LangSeedRand == 0)
    {
        LangSeedRand = rand();
        if (LangSeedRand > 1000)
        {
            LangSeedRand = LangSeedRand / 1000;
        }
    }
    if (myrank == statusrank)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 664\n");
            return 1;
        }
        fprintf(statusFile, "LangSeedRand = %d\n", LangSeedRand);
        ffcc = fclose(statusFile);
        if (ffcc != 0)
        {
            fprintf(stdout, "fclose failed 671");
            return (1);
        }
    }
    // variable commands
    sprintf(vFileNames, "variable FileNames index %s", outputfilenames);
    sprintf(vWorkPath, "variable workpath index %s", WorkPath);
    sprintf(vK, "variable K index %s", argv[6]);
    sprintf(vR0, "variable R0 index %s", argv[7]);
    sprintf(vInitFile, "variable InitFile index %s", initdataFileName);
    sprintf(vinitsteps, "variable InitSteps index %d", initdt);
    sprintf(vstepsize, "variable stepsize index %s", argv[13]);
    sprintf(vxdiff, "variable xdiff equal c_2[%d]-c_2[1]", NumBeads); // next three lines for end-to-end distance
    sprintf(vydiff, "variable ydiff equal c_3[%d]-c_3[1]", NumBeads);
    sprintf(vzdiff, "variable zdiff equal c_4[%d]-c_4[1]", NumBeads);
    sprintf(vRg, "variable rg equal c_1");
    sprintf(vLangSeed, "variable LangSeed index %d", LangSeedRand);

    // run variables in lammps
    lammps_command(lmp, vFileNames);
    lammps_command(lmp, vWorkPath);
    lammps_command(lmp, vK);
    lammps_command(lmp, vR0);
    lammps_command(lmp, vInitFile);
    lammps_command(lmp, vinitsteps);
    lammps_command(lmp, vstepsize);
    lammps_command(lmp, vxdiff);
    lammps_command(lmp, vydiff);
    lammps_command(lmp, vzdiff);
    lammps_command(lmp, vRg);
    lammps_command(lmp, vLangSeed);

    if (myrank == 0)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 708\n");
            return 1;
        }
        fprintf(statusFile, "Finished setup of variables\n");
        ffcc = fclose(statusFile);
        if (ffcc != 0)
        {
            fprintf(stdout, "fclose failed 715");
            return (1);
        }
    }

    // Run initial setup
    lammps_file(lmp, "in.FirstSteps");

    lammps_command(lmp, dumpcommand);
    lammps_command(lmp, dumpmodcommand);
    char initbondgroup[200];
    sprintf(initbondgroup, "group initbond id %d %d", 1, NumBeads); // delete old 1-N bond from previous run
    lammps_command(lmp, initbondgroup);
    lammps_command(lmp, "delete_bonds initbond bond 2 remove special");
    lammps_command(lmp, "group initbond delete");
    lammps_command(lmp, "run 0");
    lammps_command(lmp, "log none");
    MPI_Barrier(MPI_COMM_WORLD);

    // output fix. Calculates end to end distance and radius of gyration.
    char fixRendCommand[300];
    sprintf(fixRendCommand, "fix RE all ave/time 1 1 %d v_xdiff v_ydiff v_zdiff v_rg mode scalar file %sReeRg_%s.txt", dtperdump, WorkPath, outputfilenames);
    lammps_command(lmp, fixRendCommand);

    // allocate memory outside of loops for coordinates because I don't need new memory each time
    double *allcoords;
    allcoords = (double *)malloc(3 * NumBeads * sizeof(double));

    // allocate memory for contact matrices.
    int *contactmat;
    contactmat = (int *)malloc(NumBeads * NumBeads * sizeof(int));
    for (int yy = 0; yy < (NumBeads * NumBeads); yy++)
    { // initialize contactmat
        contactmat[yy] = 0;
    }

    // allocate mem for R2. one for overall calculation, one for intermittent calculation
    double *R2mat;
    double *oldR2mat;
    R2mat = (double *)malloc(NumBeads * NumBeads * sizeof(double));
    oldR2mat = (double *)malloc(sizeof(double) * NumBeads * NumBeads);
    for (int yy = 0; yy < (NumBeads * NumBeads); yy++)
    {
        R2mat[yy] = 0;
        oldR2mat[yy] = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD); // Pause everything before starting extrusion
    if (myrank == statusrank)
    {
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusfile failed line 768\n");
            return 1;
        }
        fprintf(statusFile, "Starting Extrusion...\n");
        ffcc = fclose(statusFile);
        if (ffcc != 0)
        {
            fprintf(stdout, "fclose failed 775");
            return (1);
        }
    }

    // Keep track of number of attempted moves and number of accepted moves
    int numatt[2];
    int numacc[2];
    double accperc[3];
    for (int i = 0; i < 2; i++)
    {
        numatt[i] = 0;
        numacc[i] = 0;
    }

    double addrand, UBorMove, UBrand, moveprob; // random numbers for adding, unbind or move, unbind, move
    double trytransprob = 1 - offprob;          // probability of moving vs unbinding
    int TrialLoc[2];                            // for binding

    // initiate binding number
    int bindcount;
    bindcount = -1;

    int totnumloops = 0; // keep track of which loop number i am working on
    int loopnum[maxCoh]; // which loop number
    for (int uu = 0; uu < maxCoh; uu++)
    {
        loopnum[uu] = 0;
    }

    // Start Extrusion
    for (int step = 1; step <= NumBlocks; step++)
    {                                   // for each MC step
        int numnewbonds, numdelbonds;   // initialize variables for number of new bonds and delete bonds to feed to LAMMPS
        int *newBondsAll, *delBondsAll; // arrays to hold which bonds I need to make and delete at end of MC step.

        int *newBonds; // array for new bonds. make larger just in case. Could hold zeros
        newBonds = (int *)malloc(sizeof(int) * 4 * maxCoh);
        for (int ii = 0; ii < 4 * maxCoh; ii++)
        {
            newBonds[ii] = 0;
        }
        int *delBonds; // array for deleting bonds. Could hold zeros
        delBonds = (int *)malloc(sizeof(int) * 2 * maxCoh);
        for (int ii = 0; ii < 2 * maxCoh; ii++)
        {
            delBonds[ii] = 0;
        }

        int *indwhereNew; // indices for new bons in newbonds array
        int *indwhereDel; // indices for delete bonds in delbonds array

        int yesbind[maxCoh];   // if something bound this step
        int yesunbind[maxCoh]; // if something unbinded this step
        for (int uu = 0; uu < maxCoh; uu++)
        {
            yesbind[uu] = 0;
            yesunbind[uu] = 0;
        }

        // Need to gather on all ranks
        // Get coordinates from LAMMPS. Done in run step actually, except for first step.
        if (step == 1)
        {
            lammps_gather_atoms(lmp, "x", 1, 3, allcoords);
        }

        // only do MC calculation on root rank
        if (myrank == 0)
        {
            numCoh = CurrNumCoh(cohlocs, maxCoh); // current number of bound cohesins

            int *indwherebound; // indices in the cohlifetime and cohloc arrays of the cohesins that are currently bound to chain
            indwherebound = (int *)malloc(sizeof(int) * numCoh);
            if (step == 1)
            {
                indwherebound[0] = 1;
            }
            else
            {
                indwherebound = getwherebound(indwherebound, CohLifetime, maxCoh); // current indices where cohesins are bound
            }
            int tempaddmat[NumBeads]; // temporary vector for adding cohesins. row major order
            for (int ii = 0; ii < NumBeads; ii++)
            {
                tempaddmat[ii] = 0;
            }
            int numbops = NumBeads - 1 + numCoh; // number of operations to do. NumBeads-1 = number binding attemps, numCoh=either unbind or move.
            int *whichoperation;
            whichoperation = (int *)malloc(sizeof(int) * numbops);
            whichoperation = generatepermute(whichoperation, numbops); // decide to add, move or unbind cohesins
            // 1->numCoh = unbind coh, numCoh+1->2*numCoh = move coh, 2*numCoh+1 ->
            // end = add coh
            int delIndex[maxCoh]; // 1 if deleted that coh, 0 otherwise.
            for (int iii = 0; iii < maxCoh; iii++)
            {
                delIndex[iii] = 0;
            }                        // Initialize to 0
            int indexIndelBonds = 0; // index to keep track of where to add bonds to delBonds list
            int indexInnewBonds = 0; // index to keep track of where to add bonds to newBonds list
            int AddIndex = 0;        // index to keep track of where to add bonds to tempaddmat
            int cohind;              // which cohesin I am working on;
            int currStat[2];         // bead indices to which cohesin is bound
            int thisop;              // current operation;
            int numToAdd = 0;        // number of cohesins to add this MC step

            for (int tt = 0; tt < numbops; tt++)
            {                                // for each operation
                thisop = whichoperation[tt]; // this operation. Starts from 1

                // perform operations
                if (thisop <= numCoh)
                { // Unbind or move

                    cohind = indwherebound[thisop - 1]; // which cohesin to try. starts from 1 to match LAMMPS
                    // get bead indices where bound
                    currStat[0] = cohlocs[cohind - 1];
                    currStat[1] = cohlocs[maxCoh + cohind - 1];

                    UBorMove = (double)rand() / (RAND_MAX); // decide to attempt unbinding or moving

                    // for moving
                    int currstopflag[2];
                    currstopflag[0] = stopflags[cohind - 1]; // see if it is stopped
                    currstopflag[1] = stopflags[cohind - 1 + maxCoh];
                    int *newLocs;
                    newLocs = (int *)malloc(sizeof(int) * 2);
                    // newLocs = getNewBond(newLocs, jumpprob, currStat, cohlocs, allcoords, currstopflag, NumBeads, K, R0, yesring);
                    // get new location of bond
                    newLocs = getNewBond(newLocs, jumpprob, currStat, cohlocs, allcoords, currstopflag, NumBeads, maxCoh, K, R0, yesring);
                    if (UBorMove < 0.5)
                    { // try to unbind first
                        if (currStat[0] > 0)
                        { // if something is there
                            // Unbinding code first
                            UBrand = (double)rand() / (RAND_MAX);
                            if (UBrand <= offprob)
                            { // If remove
                                // add to final life file
                                FinalLifeFile = fopen(finallifename, "a");
                                if (FinalLifeFile == NULL)
                                {
                                    fprintf(stdout, "FinalLifeFile failed line 917\n");
                                    return 1;
                                }
                                fprintf(FinalLifeFile, "%d %d\n", CohLifetime[cohind - 1], cohlocs[cohind - 1 + maxCoh] - cohlocs[cohind - 1]);
                                ffcc = fclose(FinalLifeFile);
                                if (ffcc != 0)
                                {
                                    fprintf(stdout, "fclose failed 924");
                                    return (1);
                                }
                                CohLifetime[cohind - 1] = 0; // reset coh lifetime
                                stopflags[cohind - 1] = 0;   // reset stop flags
                                stopflags[cohind - 1 + maxCoh] = 0;
                                cohlocs[cohind - 1] = 0; // reset coh locs
                                cohlocs[cohind - 1 + maxCoh] = 0;
                                delBonds[indexIndelBonds] = currStat[0];
                                delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                indexIndelBonds = indexIndelBonds + 1; // increment delBonds index
                                delIndex[cohind - 1] = 1;              // make a note that I deleted the cohesin at cohind-1;
                                yesunbind[cohind - 1] = 1;
                            }
                            else
                            { // try to move
                                if (currstopflag[0] == 0)
                                {
                                    numatt[0] = numatt[0] + 1;
                                }
                                if (currstopflag[1] == 0)
                                {
                                    numatt[1] = numatt[1] + 1;
                                }
                                if (newLocs[0] == currStat[0] && newLocs[1] == currStat[1])
                                {                                                          // if I didn't move
                                    CohLifetime[cohind - 1] = CohLifetime[cohind - 1] + 1; // increment cohlifetime
                                }
                                else if (newLocs[0] == 0)
                                { // if cohesin unbound
                                    FinalLifeFile = fopen(finallifename, "a");
                                    if (FinalLifeFile == NULL)
                                    {
                                        fprintf(stdout, "FinalLifeFile failed line 957\n");
                                        return 1;
                                    }
                                    fprintf(FinalLifeFile, "%d %d\n", CohLifetime[cohind - 1], cohlocs[cohind - 1 + maxCoh] - cohlocs[cohind - 1]);
                                    ffcc = fclose(FinalLifeFile);
                                    if (ffcc != 0)
                                    {
                                        fprintf(stdout, "fclose failed 964");
                                        return (1);
                                    }
                                    CohLifetime[cohind - 1] = 0; // reset coh lifetime
                                    stopflags[cohind - 1] = 0;   // reset stop flags
                                    stopflags[cohind - 1 + maxCoh] = 0;
                                    delIndex[cohind - 1] = 1;                // make a note that I deleted the cohesin at cohind-1;
                                    delBonds[indexIndelBonds] = currStat[0]; // populate delBonds
                                    delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                    indexIndelBonds = indexIndelBonds + 1; // increment delBonds index
                                    yesunbind[cohind - 1] = 1;
                                }
                                else
                                {                                                          // if cohesin actually moves
                                    CohLifetime[cohind - 1] = CohLifetime[cohind - 1] + 1; // increment lifetime
                                    if (indexInnewBonds == -1)
                                    {
                                        indexInnewBonds = 0;
                                    }
                                    // check if need to stop
                                    if (withstop > 0)
                                    { // if there are stop signs
                                        stopflags = updatestopflags(stopflags, currStat[0], currStat[1], newLocs, numstoppairs, leftstoplocs, rightstoplocs, cohind, maxCoh);
                                    }
                                    delBonds[indexIndelBonds] = currStat[0]; // populate delBonds for old bond
                                    delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                    indexIndelBonds = indexIndelBonds + 1;               // increment delBonds index
                                    newBonds[indexInnewBonds] = newLocs[0];              // populate newBonds
                                    newBonds[indexInnewBonds + 2 * maxCoh] = newLocs[1]; // remember that newBonds is twice as long
                                    indexInnewBonds = indexInnewBonds + 1;               // increment newBonds index;
                                                                                         // update acceptance numbers
                                    if (newLocs[0] < currStat[0])
                                    {
                                        numacc[0] = numacc[0] + 1;
                                    }
                                    if (newLocs[1] > currStat[1])
                                    {
                                        numacc[1] = numacc[1] + 1;
                                    }
                                }
                                cohlocs[cohind - 1] = newLocs[0]; // update cohlocs array
                                cohlocs[cohind - 1 + maxCoh] = newLocs[1];
                            } // end of "try to move"
                        }     // end of "if something is bound"
                    }
                    else
                    { // try to move first
                        if (currStat[0] > 0)
                        { // if something is bound
                            moveprob = (double)rand() / (RAND_MAX);
                            if (moveprob < trytransprob)
                            { // if try to move
                                if (currstopflag[0] == 0)
                                {
                                    numatt[0] = numatt[0] + 1;
                                }
                                if (currstopflag[1] == 0)
                                {
                                    numatt[1] = numatt[1] + 1;
                                }
                                if (newLocs[0] == currStat[0] && newLocs[1] == currStat[1])
                                {                                                          // if I didn't move
                                    CohLifetime[cohind - 1] = CohLifetime[cohind - 1] + 1; // increment cohlifetime
                                }
                                else if (newLocs[0] == 0)
                                { // if fell off
                                    FinalLifeFile = fopen(finallifename, "a");
                                    if (FinalLifeFile == NULL)
                                    {
                                        fprintf(stdout, "FinalLifeFile failed line 1033\n");
                                        return 1;
                                    }
                                    fprintf(FinalLifeFile, "%d %d\n", CohLifetime[cohind - 1], cohlocs[cohind - 1 + maxCoh] - cohlocs[cohind - 1]);
                                    ffcc = fclose(FinalLifeFile);
                                    if (ffcc != 0)
                                    {
                                        fprintf(stdout, "fclose failed 1040");
                                        return (1);
                                    }
                                    CohLifetime[cohind - 1] = 0; // reset coh lifetime
                                    stopflags[cohind - 1] = 0;   // reset stop flags
                                    stopflags[cohind - 1 + maxCoh] = 0;
                                    delIndex[cohind - 1] = 1;                // make a note that I deleted the cohesin at cohind-1;
                                    delBonds[indexIndelBonds] = currStat[0]; // populate delBonds
                                    delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                    indexIndelBonds = indexIndelBonds + 1; // increment delBonds index
                                    yesunbind[cohind - 1] = 1;
                                }
                                else
                                {                                                          // if I actually move
                                    CohLifetime[cohind - 1] = CohLifetime[cohind - 1] + 1; // increment lifetime
                                    if (indexInnewBonds == -1)
                                    {
                                        indexInnewBonds = 0;
                                    }
                                    // check if need to stop
                                    if (withstop > 0)
                                    { // if there are stop signs
                                        stopflags = updatestopflags(stopflags, currStat[0], currStat[1], newLocs, numstoppairs, leftstoplocs, rightstoplocs, cohind, maxCoh);
                                    }
                                    delBonds[indexIndelBonds] = currStat[0]; // populate delBonds for old bond
                                    delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                    indexIndelBonds = indexIndelBonds + 1;               // increment delBonds index
                                    newBonds[indexInnewBonds] = newLocs[0];              // populate newBonds
                                    newBonds[indexInnewBonds + 2 * maxCoh] = newLocs[1]; // remember that newBonds is twice as long
                                    indexInnewBonds = indexInnewBonds + 1;               // increment newBonds index;
                                    // update acceptance numbers
                                    if (newLocs[0] < currStat[0])
                                    {
                                        numacc[0] = numacc[0] + 1;
                                    }
                                    if (newLocs[1] > currStat[1])
                                    {
                                        numacc[1] = numacc[1] + 1;
                                    }
                                }
                                cohlocs[cohind - 1] = newLocs[0]; // update cohlocs array
                                cohlocs[cohind - 1 + maxCoh] = newLocs[1];
                            }
                            else
                            { // otherwise, unbind
                                FinalLifeFile = fopen(finallifename, "a");
                                if (FinalLifeFile == NULL)
                                {
                                    fprintf(stdout, "FinalLifeFile failed line 1088\n");
                                    return 1;
                                }
                                fprintf(FinalLifeFile, "%d %d\n", CohLifetime[cohind - 1], cohlocs[cohind - 1 + maxCoh] - cohlocs[cohind - 1]);
                                ffcc = fclose(FinalLifeFile);
                                if (ffcc != 0)
                                {
                                    fprintf(stdout, "fclose failed 1095");
                                    return (1);
                                }
                                CohLifetime[cohind - 1] = 0; // reset coh lifetime
                                stopflags[cohind - 1] = 0;   // reset stop flags
                                stopflags[cohind - 1 + maxCoh] = 0;
                                cohlocs[cohind - 1] = 0; // reset coh locs
                                cohlocs[cohind - 1 + maxCoh] = 0;
                                delBonds[indexIndelBonds] = currStat[0];
                                delBonds[indexIndelBonds + maxCoh] = currStat[1];
                                indexIndelBonds = indexIndelBonds + 1; // increment delBonds index
                                delIndex[cohind - 1] = 1;              // make a note that I deleted the cohesin at cohind-1;
                                yesunbind[cohind - 1] = 1;
                            } // end of "if move or unbind"
                        }     // end of "if something bound"
                    }         // end of "UBorMove<0.5"
                    // free memory
                    free(newLocs);
                }
                else
                { // add the cohesin
                    if (indexInnewBonds == -1)
                    {
                        indexInnewBonds = 0;
                    }
                    if (loadsite == 0)
                    {                                   // if nonspecific loading
                        int LoadBead = thisop - numCoh; // where to attempt binding. using LAMMPS indexing.
                        TrialLoc[0] = LoadBead;
                        TrialLoc[1] = LoadBead + 1; // new bond indices
                    }
                    else
                    { // specific binding
                        TrialLoc[0] = loadsite;
                        TrialLoc[1] = loadsite + 1;
                    }
                    addrand = (double)rand() / (RAND_MAX);
                    if (addrand <= onprob)
                    {                                                         // if accept addition
                        numToAdd = numToAdd + 1;                              // increment number of cohesins to add
                        tempaddmat[AddIndex] = TrialLoc[0];                   // add cohesins
                        AddIndex = AddIndex + 1;                              // increment AddIndex
                        newBonds[indexInnewBonds] = TrialLoc[0];              // populate newBonds
                        newBonds[indexInnewBonds + 2 * maxCoh] = TrialLoc[1]; // remember that newBonds is twice as long
                        indexInnewBonds = indexInnewBonds + 1;                // increment newBonds index;
                    }
                }
            } // end of "for op"
            // Finished operations for loop.

            // Now add new cohesins into cohlocs mat.
            // For now just assume that the length of all arrays is okay and does not need to change.
            if (numToAdd > 0)
            {                     // if I add cohesins
                int trackAdd = 0; // tracking which additional cohesin I am currently working on
                for (int jj = 0; jj < maxCoh; jj++)
                { // for each cohindex
                    if (CohLifetime[jj] == 0 && delIndex[jj] == 0)
                    {                                       // if nothing bound at this jj AND I did not just delete it
                        cohlocs[jj] = tempaddmat[trackAdd]; // add cohesin
                        cohlocs[jj + maxCoh] = tempaddmat[trackAdd] + 1;
                        CohLifetime[jj] = 1; // update lifetime
                        yesbind[jj] = 1;
                        loopnum[jj] = totnumloops;     // start at 0
                        totnumloops = totnumloops + 1; // keep track of which num loop overall
                        trackAdd = trackAdd + 1;       // move on to next cohesin
                        BindFile = fopen(bindfilename, "a");
                        if (BindFile == NULL)
                        {
                            fprintf(stdout, "BindFile failed line 1164\n");
                            return 1;
                        }
                        fprintf(BindFile, "%d %d %d\n", cohlocs[jj], cohlocs[jj + maxCoh], step);
                        ffcc = fclose(BindFile);
                        if (ffcc != 0)
                        {
                            fprintf(stdout, "fclose failed 1171");
                            return (1);
                        }
                        if (trackAdd == numToAdd)
                        {
                            break;
                        } // break out once I finished adding all necessary cohesins
                    }
                }
            } // end of adding

            // update number of bound cohesins and number of bonds to add/delete
            numCoh = CurrNumCoh(cohlocs, maxCoh);
            numnewbonds = CurrNumCoh(newBonds, 2 * maxCoh);
            if (numnewbonds > 0)
            {
                indwhereNew = (int *)malloc(sizeof(int) * numnewbonds);
                indwhereNew = getwherebound(indwhereNew, newBonds, 2 * maxCoh);
            }

            numdelbonds = CurrNumCoh(delBonds, maxCoh);
            if (numdelbonds > 0)
            {
                indwhereDel = (int *)malloc(sizeof(int) * numdelbonds);
                indwhereDel = getwherebound(indwhereDel, delBonds, maxCoh);
            }

            // clean up memory
            free(indwherebound);
        } // end of (if rank==0) calculations.

        // Communicate numnewbonds, numdelbonds using broadcast.
        MPI_Bcast(&numnewbonds, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numdelbonds, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&bindcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Also communicate yesbind, yesunbind, and loopnum
        MPI_Bcast(&yesbind[0], maxCoh, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&yesunbind[0], maxCoh, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&loopnum[0], maxCoh, MPI_INT, 0, MPI_COMM_WORLD);

        // setup newBondsAll array, truncate newBonds without trailing zeros at root, broadcast newBondsAll.
        // Also tell lammps to make new bonds
        if (numnewbonds > 0)
        { // if new bonds
            newBondsAll = (int *)malloc(2 * numnewbonds * sizeof(int));
            for (int oo = 0; oo < 2 * numnewbonds; oo++)
            {
                newBondsAll[oo] = 0;
            }
            // truncate without trailing zeros at root
            if (myrank == 0)
            {
                int tracknewindex = 0;
                for (int gg = 0; gg < numnewbonds; gg++)
                {
                    newBondsAll[gg] = newBonds[indwhereNew[tracknewindex] - 1];
                    newBondsAll[gg + numnewbonds] = newBonds[indwhereNew[tracknewindex] - 1 + 2 * maxCoh];
                    tracknewindex = tracknewindex + 1;
                }
                free(indwhereNew);
            }
            // broadcast newBondsAll
            MPI_Bcast(&newBondsAll[0], 2 * numnewbonds, MPI_INT, 0, MPI_COMM_WORLD);

            // Tell LAMMPS to create new bonds
            char createnewbond[100];
            for (int gg = 0; gg < numnewbonds; gg++)
            { // for each new bond
                // newbond command
                sprintf(createnewbond, "create_bonds single/bond 2 %d %d special yes", newBondsAll[gg], newBondsAll[gg + numnewbonds]);
                lammps_command(lmp, createnewbond);
            }

            // free memory
            free(newBondsAll);
        }

        // setup delBondsAll array, truncate newBonds without trailing zeros at root, broadcast delBondsAll.
        // Also tell lammps to delete bonds
        if (numdelbonds > 0)
        { // if delete bonds
            delBondsAll = (int *)malloc(2 * numdelbonds * sizeof(int));
            for (int hh = 0; hh < 2 * numdelbonds; hh++)
            {
                delBondsAll[hh] = 0;
            }
            // truncate without trailing zeros at root
            if (myrank == 0)
            {
                int trackdelindex = 0;
                for (int uu = 0; uu < numdelbonds; uu++)
                {
                    delBondsAll[uu] = delBonds[indwhereDel[trackdelindex] - 1];
                    delBondsAll[uu + numdelbonds] = delBonds[indwhereDel[trackdelindex] - 1 + maxCoh];
                    trackdelindex = trackdelindex + 1;
                }
                free(indwhereDel);
            }
            // broadcast delBondsAll;
            MPI_Bcast(&delBondsAll[0], 2 * numdelbonds, MPI_INT, 0, MPI_COMM_WORLD);

            // Tell LAMMPS to delete bonds
            char deletebondgroup[100];
            for (int uu = 0; uu < numdelbonds; uu++)
            { // for each new bond
                // delete bond commands
                sprintf(deletebondgroup, "group delbonds id %d %d", delBondsAll[uu], delBondsAll[uu + numdelbonds]);
                lammps_command(lmp, deletebondgroup);
                lammps_command(lmp, "delete_bonds delbonds bond 2 remove special");
                lammps_command(lmp, "group delbonds delete");
            }

            // free memory
            free(delBondsAll);
        }

        // Calculate run LAMMPS, gather coordinates, calculate contacts, calculate R2, calculate dist from load site to coh, calc position of "coh" COM
        // Run lammps
        lammps_command(lmp, runcommand); // runcommand defined at beginning of main. Each call = dtperMC integration steps

        // collect coords
        lammps_gather_atoms(lmp, "x", 1, 3, allcoords);

        // calculate contacts
        if (myrank == contactrank)
        {
            contactmat = getcontacts(contactmat, allcoords, NumBeads, rc);
        }
        // calculate R2
        if (myrank == R2rank)
        {
            R2mat = getR2mat(R2mat, allcoords, NumBeads);
        }

        // Update R2 matrix every contactfreq steps. Rewrite file. Rembember that we are re-averaging to make sure we don't add too much
        if (myrank == R2rank)
        {
            if ((step % contactfreq) == 0 || step == NumBlocks)
            {
                sprintf(R2name, "%sR2_%s.txt", WorkPath, outputfilenames);
                R2File = fopen(R2name, "w");
                if (R2File == NULL)
                {
                    fprintf(stdout, "R2File failed line 1314\n");
                    return 1;
                }
                for (int uu = 1; uu <= NumBeads; uu++)
                {
                    for (int nn = uu; nn <= NumBeads; nn++)
                    {
                        int thisind = (uu - 1) * NumBeads + (nn - 1);
                        double val = (oldR2mat[thisind] * (step - contactfreq) + R2mat[thisind]) / (step);
                        oldR2mat[thisind] = val;
                        R2mat[thisind] = 0; // reset R2mat
                        fprintf(R2File, "%d %d %.4f\n", uu, nn, val);
                    }
                }
                ffcc = fclose(R2File);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1331");
                    return (1);
                }
            }
        }

        // Update output contacts every contactfreq steps. Rewrite file
        if (myrank == contactrank)
        {
            if ((step % contactfreq) == 0 || step == NumBlocks)
            {
                sprintf(contactname, "%sContacts_%s.txt", WorkPath, outputfilenames);
                ContactFile = fopen(contactname, "w");
                if (ContactFile == NULL)
                {
                    fprintf(stdout, "ContactFile failed line 1346\n");
                    return 1;
                }
                for (int uu = 1; uu <= NumBeads; uu++)
                {
                    for (int nn = uu; nn <= NumBeads; nn++)
                    {
                        int thiselement = contactmat[(uu - 1) * NumBeads + (nn - 1)];
                        fprintf(ContactFile, "%d %d %d\n", uu, nn, thiselement);
                    }
                }
                ffcc = fclose(ContactFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1360\n");
                    return (1);
                }
            }
        }

        // Print status updates from rank 1 (if exists), or 0 (defined by statusrank variable at begining of main)
        if (myrank == statusrank)
        {
            if ((step % 10000) == 0)
            {
                statusFile = fopen(statusfilename, "a");
                if (statusFile == NULL)
                {
                    fprintf(stdout, "statusFile failed line 1374\n");
                    return 1;
                }
                fprintf(statusFile, "Block %d\n", step);
                ffcc = fclose(statusFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1381\n");
                    return (1);
                }
            }
            if (step == round(0.25 * NumBlocks))
            {
                statusFile = fopen(statusfilename, "a");
                if (statusFile == NULL)
                {
                    fprintf(stdout, "statusFile failed line 1390\n");
                    return 1;
                }
                fprintf(statusFile, "**** 25 percent done ****\n");
                ffcc = fclose(statusFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1397");
                    return (1);
                }
            }
            else if (step == round(0.5 * NumBlocks))
            {
                statusFile = fopen(statusfilename, "a");
                if (statusFile == NULL)
                {
                    fprintf(stdout, "statusFile failed line 1406\n");
                    return 1;
                }
                fprintf(statusFile, "**** 50 percent done ****\n");
                ffcc = fclose(statusFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1413");
                    return (1);
                }
            }
            else if (step == round(0.75 * NumBlocks))
            {
                statusFile = fopen(statusfilename, "a");
                if (statusFile == NULL)
                {
                    fprintf(stdout, "statusFile failed line 1422\n");
                    return 1;
                }
                fprintf(statusFile, "**** 75 percent done ****\n");
                ffcc = fclose(statusFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1429");
                    return (1);
                }
            }
            else if (step == (NumBlocks - 1))
            {
                statusFile = fopen(statusfilename, "a");
                if (statusFile == NULL)
                {
                    fprintf(stdout, "statusFile failed line 1438\n");
                    return 1;
                }
                fprintf(statusFile, "**** Last step ****\n");
                ffcc = fclose(statusFile);
                if (ffcc != 0)
                {
                    fprintf(stdout, "fclose failed 1445");
                    return (1);
                }
            }
        }

        // print number of cohesins, number of attempted moves, number of accepted moves from root
        if (myrank == 0)
        {

            // numcoh
            NumCohFile = fopen(numcohname, "a");
            if (NumCohFile == NULL)
            {
                fprintf(stdout, "NumCohFile failed line 1459\n");
                return 1;
            }
            fprintf(NumCohFile, "%d\n", numCoh);
            ffcc = fclose(NumCohFile);
            if (ffcc != 0)
            {
                fprintf(stdout, "fclose failed 1466");
                return (1);
            }

            // Numatt and Numacc
            AccFile = fopen(accfilename, "a");
            if (AccFile == NULL)
            {
                fprintf(stdout, "AccFile failed line 1474\n");
                return 1;
            }
            fprintf(AccFile, "%d %d %d %d\n", numatt[0], numacc[0], numatt[1], numacc[1]);
            ffcc = fclose(AccFile);
            if (ffcc != 0)
            {
                fprintf(stdout, "fclose failed 1481");
                return (1);
            }
        } // Finished printing

        free(newBonds);
        free(delBonds);
        MPI_Barrier(MPI_COMM_WORLD); // hold at end of single MC step
    }                                // end of MC steps
    if (myrank == statusrank)
    { // print percent accepted
        accperc[0] = ((double)numacc[0]) / ((double)numatt[0]);
        accperc[1] = ((double)numacc[1]) / ((double)numatt[1]);
        accperc[2] = ((double)numacc[0] + (double)numacc[1]) / ((double)numatt[0] + (double)numatt[1]);
        statusFile = fopen(statusfilename, "a");
        if (statusFile == NULL)
        {
            fprintf(stdout, "statusFile failed line 1498\n");
            return 1;
        }
        fprintf(statusFile, "Percent accepted = \n Left domain: %.3f\n Right domain: %.3f\n Overall: %.3f \n", accperc[0], accperc[1], accperc[2]);
        ffcc = fclose(statusFile);
        if (ffcc != 0)
        {
            fprintf(stdout, "fclose failed 1505");
            return (1);
        }
    }

    lammps_close(lmp);
    free(allcoords);
    free(cohlocs);
    free(CohLifetime);
    free(leftstoplocs);
    free(rightstoplocs);
    free(stopflags);
    free(R2mat);
    free(oldR2mat);
    free(contactmat);
    // free(CohCOM);
    MPI_Finalize();
    return 0;
} // end of main
