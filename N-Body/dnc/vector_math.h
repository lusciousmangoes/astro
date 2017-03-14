
#ifndef VECTOR_CLEAR
#define VECTOR_CLEAR(v){(v)[0] = 0.0; (v)[1] = 0.0; (v)[2] = 0.0;}
#endif

#ifndef VECTOR6_CLEAR
#define VECTOR6_CLEAR(v){(v)[0] = 0.0; (v)[1] = 0.0; (v)[2] = 0.0;\
                         (v)[3] = 0.0; (v)[4] = 0.0; (v)[5] = 0.0;}
#endif

#ifndef VECTOR10_CLEAR
#define VECTOR10_CLEAR(v){(v)[0] = 0.0; (v)[1] = 0.0; (v)[2] = 0.0;\
                          (v)[3] = 0.0; (v)[4] = 0.0; (v)[5] = 0.0;\
                          (v)[6] = 0.0; (v)[7] = 0.0; (v)[8] = 0.0;\
                          (v)[9] = 0.0; }
#endif

#ifndef VECTOR_SUB
#define VECTOR_SUB(d, v1, v2){ \
                    (d)[0]=(v2)[0]-(v1)[0]; \
                    (d)[1]=(v2)[1]-(v1)[1]; \
                    (d)[2]=(v2)[2]-(v1)[2]; }
#endif

#ifndef VECTOR10_INCREMENT
#define VECTOR10_INCREMENT(d, v){ \
                    (d)[0]+=(v)[0]; \
                    (d)[1]+=(v)[1]; \
                    (d)[2]+=(v)[2]; \
                    (d)[3]+=(v)[3]; \
                    (d)[4]+=(v)[4]; \
                    (d)[5]+=(v)[5]; \
                    (d)[6]+=(v)[6]; \
                    (d)[7]+=(v)[7]; \
                    (d)[8]+=(v)[8]; \
                    (d)[9]+=(v)[9]; }
#endif

#ifndef VECTOR_SET
#define VECTOR_SET(b,a){(b)[0] = (a)[0]; (b)[1] = (a)[1]; (b)[2] = (a)[2];}
#endif

#ifndef VECTOR_MAG2
#define VECTOR_MAG2(m2,v){(m2) = (v)[0]*(v)[0] + (v)[1]*(v)[1] + (v)[2]*(v)[2];}
#endif

#ifndef VECTOR_MUL_ADD
#define VECTOR_MUL_ADD(b, m, a){ \
                    (b)[0] += (m) * (a)[0]; \
                    (b)[1] += (m) * (a)[1]; \
                    (b)[2] += (m) * (a)[2]; }
#endif

#ifndef VECTOR6_MUL_ADD
#define VECTOR6_MUL_ADD(b, m, a){ \
                    (b)[0] += (m) * (a)[0]; \
                    (b)[1] += (m) * (a)[1]; \
                    (b)[2] += (m) * (a)[2]; \
                    (b)[3] += (m) * (a)[3]; \
                    (b)[4] += (m) * (a)[4]; \
                    (b)[5] += (m) * (a)[5]; }
#endif

#ifndef VECTOR10_MUL_ADD
#define VECTOR10_MUL_ADD(b, m, a){ \
                    (b)[0] += (m) * (a)[0]; \
                    (b)[1] += (m) * (a)[1]; \
                    (b)[2] += (m) * (a)[2]; \
                    (b)[3] += (m) * (a)[3]; \
                    (b)[4] += (m) * (a)[4]; \
                    (b)[5] += (m) * (a)[5]; \
                    (b)[6] += (m) * (a)[6]; \
                    (b)[7] += (m) * (a)[7]; \
                    (b)[8] += (m) * (a)[8]; \
                    (b)[9] += (m) * (a)[9]; }
#endif

#ifndef VECTOR_MUL
#define VECTOR_MUL(v, m){(v)[0]*=(m); (v)[1]*=(m);(v)[2]*=(m);}
#endif
#ifndef VECTOR6_MUL
#define VECTOR6_MUL(v, m){(v)[0]*=(m); (v)[1]*=(m);(v)[2]*=(m); \
                          (v)[3]*=(m); (v)[4]*=(m);(v)[5]*=(m);}
#endif

#ifndef VECTOR_DIV
#define VECTOR_DIV(v, d){ float _recip; _recip=1.0/(d); VECTOR_MUL((v), _recip); }
#endif

#ifndef VECTOR6_DIV
#define VECTOR6_DIV(v, d){ float _recip; _recip=1.0/(d); VECTOR6_MUL((v), _recip); }
#endif

#ifndef RECIP_ROOT
#define RECIP_ROOT(rr,d2) (rr) = 1.0/sqrt((d2))
#endif

