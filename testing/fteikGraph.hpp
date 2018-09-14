#ifndef FTEIK_GRAPH_HPP__
#define FTEIK_GRAPH_HPP__ 1

class fteikGraph
{
    public:
       fteikGraph(const int nz, const int nx, const int ny, int verbose=0);
       ~fteikGraph(void);
       size_t getNumberOfGridPoints(void);
       void finalize(void);
       int computeNodeLevels(const int sweep);
       int *getLevelPointer(const int sweep) __attribute__((aligned(64)));
       int *getLevelPtrPointer(const int sweep) __attribute__((aligned(64)));
       int copyIJKV(const int nwork, const int sweep,
                    int *__restrict__ ijkvOut);
       int *getIJKVPointer(const int sweep) __attribute__((aligned(64)));
       bool isInitialized(void);
       int copyLevelPtr(const int nwork, const int sweep,
                        int *__restrict__ levelPtrOut);
/*
       int generateLevelpointers(const size_t ngrd, const int nLevels,
                                 const int *__restrict__ levels,
                                 int *__restrict__ levelPtr,
                                 int *__restrict__ ijkv);
*/
       void sweep2sgnx64f(const int sweep,
                          int *sgntz, int *sgntx, int *sgnty,
                          int *sgnvz, int *sgnvx, int *sgnvy,
                          double *sgnrz, double *sgnrx, double *sgnry);
       int grid2index(const int i, const int j, const int k)
       {
           return k*nz*nx + j*nz + i;
       };
       int index2grid(const int igrd,
                      int *i, int *j, int *k)
       {
           int ierr, nzx;
           ierr = 0;
           nzx = nz*nx;
           *k = (igrd)/nzx;
           *j = (igrd - *k*nzx)/nz;
           *i =  igrd - *k*nzx - *j*nz;
           if (*i < 0 || *i > nz - 1){ierr = ierr + 1;}
           if (*j < 0 || *j > nx - 1){ierr = ierr + 1;}
           if (*k < 0 || *k > ny - 1){ierr = ierr + 1;}
           return ierr;
       };
       int getMaxLevelSize(void)
       {
           return maxLevelSize;
       };
       int getNumberOfLevels(void);

    private:
       int generateIJKV(const int sweep);
       void allocateLevelPtrs(void);
       void setSweepSigns(void);
       void generateSweepOrderings(void);
       int *ijkv1, *ijkv2, *ijkv3, *ijkv4,
           *ijkv5, *ijkv6, *ijkv7, *ijkv8;
       int *levels1;
/*
       int *levels2, *levels3, *levels4,
           *levels5, *levels6, *levels7, *levels8;
*/
       int *levelPtr1;
/*
       int *levelPtr2, *levelPtr3, *levelPtr4,
           *levelPtr5, *levelPtr6, *levelPtr7, *levelPtr8;
*/
       int nz, nx, ny;
       int nLevels;
       int maxLevelSize;
       int sgntzSweep1, sgntxSweep1, sgntySweep1,
           sgntzSweep2, sgntxSweep2, sgntySweep2,
           sgntzSweep3, sgntxSweep3, sgntySweep3,
           sgntzSweep4, sgntxSweep4, sgntySweep4,
           sgntzSweep5, sgntxSweep5, sgntySweep5,
           sgntzSweep6, sgntxSweep6, sgntySweep6,
           sgntzSweep7, sgntxSweep7, sgntySweep7,
           sgntzSweep8, sgntxSweep8, sgntySweep8;
       int sgnvzSweep1, sgnvxSweep1, sgnvySweep1,
           sgnvzSweep2, sgnvxSweep2, sgnvySweep2,
           sgnvzSweep3, sgnvxSweep3, sgnvySweep3,
           sgnvzSweep4, sgnvxSweep4, sgnvySweep4,
           sgnvzSweep5, sgnvxSweep5, sgnvySweep5,
           sgnvzSweep6, sgnvxSweep6, sgnvySweep6,
           sgnvzSweep7, sgnvxSweep7, sgnvySweep7,
           sgnvzSweep8, sgnvxSweep8, sgnvySweep8;
       double sgnrzSweep1, sgnrxSweep1, sgnrySweep1,
              sgnrzSweep2, sgnrxSweep2, sgnrySweep2,
              sgnrzSweep3, sgnrxSweep3, sgnrySweep3,
              sgnrzSweep4, sgnrxSweep4, sgnrySweep4,
              sgnrzSweep5, sgnrxSweep5, sgnrySweep5,
              sgnrzSweep6, sgnrxSweep6, sgnrySweep6,
              sgnrzSweep7, sgnrxSweep7, sgnrySweep7,
              sgnrzSweep8, sgnrxSweep8, sgnrySweep8; 
       int z1Sweep1, z1Sweep2, z1Sweep3, z1Sweep4,
           z1Sweep5, z1Sweep6, z1Sweep7, z1Sweep8,
           z2Sweep1, z2Sweep2, z2Sweep3, z2Sweep4,
           z2Sweep5, z2Sweep6, z2Sweep7, z2Sweep8;
       int x1Sweep1, x1Sweep2, x1Sweep3, x1Sweep4,
           x1Sweep5, x1Sweep6, x1Sweep7, x1Sweep8,
           x2Sweep1, x2Sweep2, x2Sweep3, x2Sweep4,
           x2Sweep5, x2Sweep6, x2Sweep7, x2Sweep8;
       int y1Sweep1, y1Sweep2, y1Sweep3, y1Sweep4,
           y1Sweep5, y1Sweep6, y1Sweep7, y1Sweep8,
           y2Sweep1, y2Sweep2, y2Sweep3, y2Sweep4,
           y2Sweep5, y2Sweep6, y2Sweep7, y2Sweep8;
       size_t ngrd; /*!< Number of grid points (=nz*nx*ny) */
       size_t nnz;  /*!< Number of non-zeros in global matrix. */
       int verbose;
       bool linit;
};

#endif
