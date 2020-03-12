#ifndef __MT19937AR_INCLUDED__
#define __MT19937AR_INCLUDED__

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,1)-real-interval */
double genrand(void);

#endif
