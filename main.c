#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define MAX(a,b) (a > b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.iv]
#define RECURSIVE 0

typedef struct node {
	int k;
	int ix, iy, iv;
	double x, y, v;
} node_t;

typedef struct {
	double u;
	int dy;
} ctrl_t;

double ****J;
ctrl_t ****C;

/* problem */
struct {
	double vd, T, C;
	double x0, y0, v0;
	double *obst_x, *obst_y, *obst_v;
	int obst_x_len, obst_y_len, obst_v_len;
	int numlanes, numsteps, n;
}P = {0};

/* solution */
struct {
	double **obst_x;
	double **obst_y;
	double **obst_v;
}S = {0};

/* state variable domains */
struct {
	int NX;
	double dx, *X;

	int NY;
	double *Y;

	int NV;
	double vmin, vmax, dv, *V;

	double umin, umax;
}D = {0};

static void
read_P_addval(char *buf, char *name, double **dst, int *dstlen)
{
	char fmt[1024];
	int ival;
	double dval;

	snprintf(fmt, sizeof(fmt), "\"%s(%s)\":%s\n", name, "%d", "%lf");
	if (sscanf(buf, fmt, &ival, &dval)) {
		if (*dstlen <= ival) {
			*dstlen = ival;
			assert(((*dst) = realloc((*dst), sizeof(double)*(*dstlen))));
		}
		(*dst)[ival-1] = dval;
	} 
}

/* load problem in structure P */
static void
read_P(void)
{
	char buf[1024];
	double dval;
	int ival;
	while(fgets(buf, sizeof(buf), stdin))
	{
		if (sscanf(buf, "\"vd\":%lf\n", &dval) == 1) 
			P.vd = dval;
		if (sscanf(buf, "\"numlanes\":%d\n", &ival) == 1) 
			P.numlanes = ival;
		if (sscanf(buf, "\"numsteps\":%d\n", &ival) == 1) 
			P.numsteps = ival;
		if (sscanf(buf, "\"C\":%lf\n", &dval) == 1) 
			P.C = dval;
		if (sscanf(buf, "\"T\":%lf\n", &dval) == 1) 
			P.T = dval;
		if (sscanf(buf, "\"x(0)\":%lf\n", &dval) == 1) 
			P.x0 = dval;
		if (sscanf(buf, "\"y(0)\":%lf\n", &dval) == 1) 
			P.y0 = dval;
		if (sscanf(buf, "\"v(0)\":%lf\n", &dval) == 1) 
			P.v0 = dval;
		read_P_addval(buf, "obst_x", &P.obst_x, &P.obst_x_len);
		read_P_addval(buf, "obst_y", &P.obst_y, &P.obst_y_len);
		read_P_addval(buf, "obst_v", &P.obst_v, &P.obst_v_len);
	}

	if (!(P.obst_x_len == P.obst_y_len && P.obst_x_len == P.obst_v_len)) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	P.n = P.obst_x_len;
}

/* compute obstacle trajectories in structure S */
static void
init_S(void)
{
	int i, k;

	/* allocate */
	assert((S.obst_x = calloc(sizeof(double*), P.numsteps)));
	assert((S.obst_y = calloc(sizeof(double*), P.numsteps)));
	assert((S.obst_v = calloc(sizeof(double*), P.numsteps)));
	for (k=0; k < P.numsteps; k++) {
		assert((S.obst_x[k] = calloc(sizeof(double), P.n)));
		assert((S.obst_y[k] = calloc(sizeof(double), P.n)));
		assert((S.obst_v[k] = calloc(sizeof(double), P.n)));
	}

	/* compute */
	for (k=0; k < P.numsteps; k++) {
		if (k==0) {
			for (i=0; i < P.n; i++) {
				S.obst_x[k][i] = P.obst_x[i];
				S.obst_y[k][i] = P.obst_y[i];
				S.obst_v[k][i] = P.obst_v[i];
			}
			continue;
		}
		for (i=0; i < P.n; i++) {
			S.obst_x[k][i] = S.obst_x[k-1][i] + S.obst_v[k-1][i]*P.T;
			S.obst_y[k][i] = S.obst_y[k-1][i];
			S.obst_v[k][i] = S.obst_v[k-1][i];
		}
	}
}

static void
init_D(void)
{
	int ix, iy, iv;

	/* x domain */
	D.NX = 1024;
	D.dx = (1.5*P.vd * P.T * P.numsteps)/D.NX;
	if (D.dx >= 1.0*P.C) {
		D.dx = 0.5*P.C;
		D.NX = (int)ceil((1.5*P.vd * P.T * P.numsteps)/D.dx);
	}
	assert((D.X = calloc(sizeof(double), D.NX)));
	for (ix = 0; ix < D.NX; ix++) {
		D.X[ix] = P.x0 + D.dx * ix;
		fprintf(stderr, "%.1f ", D.X[ix]);
	}
	fprintf(stderr, "\n");

	/* y domain */
	D.NY = P.numlanes;
	assert((D.Y = calloc(sizeof(double), D.NY)));
	for (iy = 0; iy < D.NY; iy++)
		D.Y[iy] = iy;

	/* v domain */
	D.NV = 20;
	D.vmin = 0;
	D.vmax = MAX(1.1*P.v0, 2*P.vd);
	D.dv = (D.vmax - D.vmin)/D.NV;
	assert((D.V = calloc(sizeof(double), D.NV)));
	for (iv = 0; iv < D.NV; iv++) 
		D.V[iv] = D.dv * iv;

	D.umin = -12;
	D.umax = 8;
}

static void 
init_JC(void)
{
	assert((J = calloc(sizeof(double*), P.numsteps)));
	assert((C = calloc(sizeof(ctrl_t*), P.numsteps)));
	int k;
	for (k=0; k < P.numsteps; k++) {
		assert((J[k] = calloc(sizeof(double*), D.NX)));
		assert((C[k] = calloc(sizeof(ctrl_t*), D.NX)));
		int ix;
		for (ix=0; ix < D.NX; ix++) {
			assert((J[k][ix] = calloc(sizeof(double*), D.NY)));
			assert((C[k][ix] = calloc(sizeof(ctrl_t*), D.NY)));
			int iy;
			for (iy=0; iy < D.NY; iy++) {
				assert((J[k][ix][iy] = calloc(sizeof(double), D.NV)));
				assert((C[k][ix][iy] = calloc(sizeof(ctrl_t), D.NV)));
			}
		}
	}
}

static int
discretizexv(node_t *state)
{
	int ix, iv;

	ix = (int)round((state->x - D.X[0])/D.dx);
	iv = (int)round((state->v - D.V[0])/D.dv);

	if (!(ix >= 0 && ix < D.NX))
		return 1;
	if (!(iv >= 0 && iv < D.NV))
		return 2;

	state->x = D.X[(state->ix = ix)];
	state->v = D.V[(state->iv = iv)];

	return 0;
}

static double
crash(node_t state)
{
	int i;
	for (i=0; i < P.n; i++) {
		if (S.obst_y[state.k][i] != state.y)
			continue;
		if (fabs(S.obst_x[state.k][i] - state.x) <= 1.0*P.C)
			return 10e5;
	}
	return 0;
}

#if(RECURSIVE)
static void
recursive_dp(node_t this)
{
	node_t next = {0};
	double dy, u, phi;
	double *Jthis, *Jnext;
	ctrl_t *Cthis;

	if (this.k == 0)
		assert(this.ix == 0 && this.iy == 0 && this.iv == 0);

	/* final stage costs */
	if (this.k == P.numsteps-1) {
		J[this.k][this.ix][this.iy][this.iv] = crash(this) + pow(this.v - P.vd, 2);
		return;
	}

	/* intermediate state */
	Jthis = &J[this.k][this.ix][this.iy][this.iv];
	Cthis = &C[this.k][this.ix][this.iy][this.iv];

	*Jthis = DBL_MAX;
	next.k = this.k+1;
	for (u=D.umin; u <= D.umax; u+=4) {
		next.x = this.x + this.v*P.T + 0.5*u*pow(P.T, 2);
		next.v = this.v + u*P.T;

		if (discretizexv(&next))
			continue;

		for (dy=-1; dy <= 1; dy++) {
			next.y = this.y + dy;
			if (!(next.y >= 0 && next.y < D.NY))
				continue;
			next.iy = (int)next.y;
			fprintf(stderr, "next.y %.1f next.iy %d\n", next.y, next.iy);

			Jnext = &J[next.k][next.ix][next.iy][next.iv];

			if (*Jnext == 0)
				recursive_dp(next);


			phi = crash(this) + pow(u, 2) + pow(dy, 2) + pow(next.v - P.vd, 2);

			if (phi + *Jnext < *Jthis) {
				*Jthis = phi + *Jnext;
				Cthis->u = u;
				Cthis->dy = dy;
			}
		}
	}

	//printf("%d %d %d %d %.1f\n", this.k, this.ix, this.iy, this.iv, *Jthis);
}
#endif

static node_t *
getnext(node_t this, ctrl_t c)
{
	static node_t next;

	next.k = this.k+1;
	next.x = this.x + this.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = this.v + c.u*P.T;
	if (discretizexv(&next)) {
		fprintf(stderr, "getnext: discretize failed\n");
		return NULL;
	}

	next.y = this.y + c.dy;
	next.iy = (int)next.y;
	if (next.y < 0 || next.y >= D.NY) {
		return NULL;
	}


	return &next;
}

static void
iterative_dp(void)
{
	int k, ix, iy, iv;
	double x, y, v;
	for (k=P.numsteps-1; k >= 0; k--) {
		for (ix=0; ix < ((k == 0)?1:D.NX); ix++) {
			x = ((k == 0)?P.x0:D.X[ix]);
			for (iy=0; iy < ((k == 0)?1:D.NY); iy++) {
				y = ((k == 0)?P.y0:D.Y[iy]);
				for (iv=0; iv < ((k == 0)?1:D.NV); iv++) {
					v = ((k == 0)?P.v0:D.V[iv]);

					node_t this = {
						.k = k, 
						.x = x, 
						.y = y, 
						.v = v, 
						.ix = ((k == 0)?0: ix), 
						.iy = ((k == 0)?0: (int)y), 
						.iv = ((k == 0)?0: iv)
					};

					if (this.ix > D.NX || this.iy > D.NY || this.iv > D.NV) {
						fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f y %.1f/%.1f v %.1f/%.1f\n",
								this.k, this.x, D.X[D.NX-1], this.y, 
								(double)P.numlanes-1, this.v, D.V[D.NV-1]);
						exit(1);
					}

					node_t next = {0};
					double dy, u, phi;
					double *Jthis, *Jnext;
					ctrl_t *Cthis;

					if (this.k == 0) 
						assert(this.ix == 0 && this.iy == 0 && this.iv == 0);

					if (this.k == P.numsteps-1) {
						/* final stage costs */
						J[this.k][this.ix][this.iy][this.iv] = crash(this) + pow(this.v - P.vd, 2);
						continue;
					}

					/* intermediate state */
					Jthis = &J[this.k][this.ix][this.iy][this.iv];
					Cthis = &C[this.k][this.ix][this.iy][this.iv];

					*Jthis = DBL_MAX;
					next.k = this.k+1;
					for (u=D.umin; u <= D.umax; u+=4) {
						next.x = this.x + this.v*P.T + 0.5*u*pow(P.T, 2);
						next.v = this.v + u*P.T;

						if (discretizexv(&next))
							continue;

						for (dy=-1; dy <= 1; dy++) {
							next.y = this.y + dy;
							if (!(next.y >= 0 && next.y < D.NY))
								continue;

							next.iy = (int)next.y;
							assert(next.y == (double)next.iy);

							Jnext = &J[next.k][next.ix][next.iy][next.iv];

							phi = crash(this) + pow(this.v - P.vd, 2) + pow(u, 2) + pow(dy, 2);

							if (phi + *Jnext < *Jthis) {
								*Jthis 		= phi + *Jnext;
								Cthis->u		= u;
								Cthis->dy	= dy;
							}
						}
					}
				}
			}
		}
	}
}

static void
printsol(node_t initial)
{
	node_t *nextptr;
	int k = 0, i;
	node_t state[P.numsteps];


	printf("// cost = %.1f\n", J[0][0][0][0]);
	puts("var Data = {");

	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c;
	for (k=0; k < P.numsteps-1; k++) {
		c = ACCESS(C, state[k]);
		fprintf(stderr, "k %d state %.1f(%d) %.1f(%d) %.1f(%d) \t J %.1f \t C %.1f %d -->\n", 
				k, state[k].x, state[k].ix, 
				state[k].y, state[k].iy, 
				state[k].v, state[k].iv, 
				ACCESS(J, state[k]),
				c.u, c.dy);
		if ((nextptr = getnext(state[k], c)) == NULL)  {
			fprintf(stderr, "invalid \n");
			exit(1);
		}
		state[k+1] = *nextptr;
	} 

	/* ego trajectory from squence of states/controls */
	double x = state[0].x;
	double v = state[0].v;
	double u, dy;
	for (k=0; k < P.numsteps; k++) {
		u = ACCESS(C, state[k]).u;
		dy = ACCESS(C, state[k]).dy;
		printf("\"x(%d)\":%.1f, /*%.1f*/\n", k, state[k].x, state[k].x);
		printf("\"y(%d)\":%.1f,\n", k, state[k].y);
		printf("\"vx(%d)\":%.1f, /*%.1f*/\n", k, v, state[k].v);
		printf("\"vy(%d)\":%.1f,\n", k, 0.0);
		printf("\"ux(%d)\":%.1f, ", k, u);
		printf("\"uy(%d)\":%.1f,\n", k, dy);
		puts("");
		x = x + (v * P.T) + (0.5 * u * pow(P.T, 2));
		v = v + (u * P.T);
	}

	/* obstacle trajectories */
	for (k=0; k < P.numsteps; k++)  {
		for (i=0; i < P.n; i++) {
			printf("\"obst_x(%d,%d)\":%.1f, ", i+1, k, S.obst_x[k][i]);
			printf("\"obst_y(%d,%d)\":%.1f,\n", i+1, k, S.obst_y[k][i]);
		}
		puts("");
	}
	
	printf("\"n\":%d,", P.n);
	printf("\"k\":%d,", P.numsteps);
	printf("\"numlanes\":%d,", P.numlanes);
	printf("\"vd\":%f,", P.vd);
	printf("\"Step\":%f", P.T);
	printf("};\n");
}

int
main(void)
{
	read_P();
	init_S();
	init_D();
	init_JC();

	fprintf(stderr, "P.numlanes: %d %d\n", P.numlanes, D.NY);

	node_t initial = {
		.k = 0,
		.x = P.x0,
		.y = P.y0,
		.v = P.v0
	};

#if(RECURSIVE)
	recursive_dp(initial);
#else
	iterative_dp();
#endif

	printsol(initial);

	return 0;
}

