
#if WORDSIZE == 32
#define lsode dlsode
#define b22in db22in
#define bs2dr dbs2dr
#define bs2vl dbs2vl
#define bs2in dbs2in
#define bs3in dbs3in
#define bs3vl dbs3vl
#define bsnak dbsnak
#define isamax idamax_u
#define sasum dasum_u
#define saxpy daxpy_u
#define scopy dcopy_u
#define sdot ddot_u
#define sgbco dgbco_u
#define sgbfa dgbfa_u
#define sgbsl dgbsl_u
#define sgefa dgefa_u
#define sgesl dgesl_u
#define snrm2 dnrm2_u
#define sscal dscal_u
#define ssum dsum
#define sswap dswap_u
#define ISAMAX idamax_u
#define SASUM dasum_u
#define SAXPY daxpy_u
#define SCOPY dcopy_u
#define SDOT ddot_u
#define SSCAL dscal_u
#define SSWAP dswap_u
#define SNRM2 dnrm2_u
#define SSUM dsum
#define SGBFA dgbfa_u
#define SGBSL dgbsl_u
#define SGEFA dgefa_u
#define SGESL dgesl_u
#define EWSET dewset_u
#define XERRWV dxerrwv_u
#define dfloat real
#define float real
#endif