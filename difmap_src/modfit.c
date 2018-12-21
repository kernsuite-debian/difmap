#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "obs.h"
#include "vlbconst.h"
#include "logio.h"
#include "lmfit.h"
#include "besj.h"

typedef struct { /* The partial derivative of the model vs one free parameter */
  double re;     /* Real part of partial derivative */
  double im;     /* Imaginary part of partial derivative */
} Vispar;

/*
 * Define an object to contain state info during fitting.
 */
typedef struct {
  Lmfit *lm;       /* The Levensberg-Marquardt fit object */
  Observation *ob; /* The descriptor of the observation */
  Model *mod;      /* The model being fitted */
  int nfree;       /* The number of free parameters in 'mod' */
  int cif;         /* Index of the IF to be processed next */
  int isub;        /* Index of the sub-array to be processed next */
  int itime;       /* Index of the integration from isub to be processed next */
  int ibase;       /* Index of the baseline from itime to be processed next */
  int done;        /* Start new visibility next? */
  int eod;         /* True when no more data require processing */
  Vispar *vp;      /* Array of nfree visibility complex partial derivatives */
  double re;       /* Real part of (data - model) */
  double im;       /* Imaginary part of (data - model) */
  double wt;       /* Weight of visibility */
  double uvrmin;   /* The minimum UV radius */
  double uvrmax;   /* The maximum UV radius - used to renormalize U and V */
} Modfit;

static Modfit *new_Modfit(Observation *ob, Model *mod,
			  float uvmin, float uvmax);
static Modfit *del_Modfit(Modfit *mf);
static int mod_nfree(Model *mod);
static int endfit(Observation *ob, Modfit *mf, int quiet, int retcode);
static int getmodvis(Modfit *mf, Visibility *vis);
static int update_model_errors(Modfit *mf);

static GETFREE(getfree);
static SETFREE(setfree);
static GETNEXT(getnext);

/*.......................................................................
 * Fit the variable components of the established and tentative models to
 * the visibilities of the currently selected processing stream,
 * via non-linear least-squares fitting to the real and imaginary
 * parts of the residual visibilities in the UV plane.
 *
 * Input:
 *  ob   Observation *  The observation to fit the model to.
 *  niter        int    The number of iterations to perform, or 0 to
 *                      request just the pre-pass iteration.
 *  uvmin      float    The minimum UV radius to take visibilities from.
 *  uvmax      float    The maximum UV radius to take visibilities from.
 *  quiet        int    0 - Display the results of the fit to standard output.
 *                      1 - Don't display the results of the fit.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int fituvmodel(Observation *ob, int niter, float uvmin, float uvmax, int quiet)
{
  Modfit *mf;     /* Model-fit descriptor */
  Lmfit *lm;      /* Levensburgh-Marquardt fit descriptor */
  int iter;       /* The iteration being performed */
  int old_if;     /* The current IF state to be restored on exit */
  int was_best=1; /* True when the previous iteration was the best fit so far */
/*
 * Valid observation?
 */
  if(!ob_ready(ob, OB_SELECT, "fituvmodel"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Move all fixed components into the established model, and place the
 * variable components in the tentative model.
 */
  if(obvarmod(ob))
    return 1;
/*
 * Create the model-fit descriptor.
 */
  mf = new_Modfit(ob, ob->newmod, uvmin, uvmax);
  if(mf==NULL)
    return endfit(ob, mf, quiet, 1);
/*
 * Get the Levensberg-Marquardt fit descriptor.
 */
  lm = mf->lm;
/*
 * Fit the model.
 */
  for(iter=0; iter<=niter; iter++) {
/*
 * Perform another iteration.
 */
    switch(lm_fit(lm)) {
    case LM_ABORT:
      return endfit(ob, mf, quiet, 1);
      break;
    case LM_BETTER:
/*
 * If the previous lines were a block of failed iterations, add a blank line
 * to end the block.
 */
      if(!was_best && !quiet)
	lprintf(stdout, "\n");
/*
 * Describe the scope of the fitting problem on the first iteration.
 */
      if(iter==0 && !quiet) {
	long nvis = (lm->best.ndfree + lm->nfree) / 2L;
	lprintf(stdout,"There are %d variables and %ld usable visibilities.\n",
		lm->nfree, nvis);
	lprintf(stdout, "This gives 2 x %ld - %d = %ld degrees of freedom.\n",
		nvis, lm->nfree, lm->best.ndfree);
	lprintf(stdout, "Reduced Chi-squared = Chi-squared / %ld.\n\n",
		lm->best.ndfree);
      };
/*
 * Report the improved fit.
 */
      if(!quiet) {
        lprintf(stdout,
                "Iteration %2.2d: Reduced Chi-squared=%#.8g  Degrees of Freedom=%ld\n",
                iter, lm->best.rchisq, lm->best.ndfree);
        wmodel(ob->newmod, 0.0f, 0.0f, 0, 0.0f, stdout);

/*
 * Separate the fit details from those of the next iteration.
 */
        lprintf(stdout, "\n");
      }
/*
 * Record the fact that a new best fit was attained.
 */
      was_best = 1;
      break;
    case LM_WORSE:
/*
 * The latest iteration did not preduce a better fit.
 */
      if(!quiet) {
        lprintf(stdout,"Iteration %2.2d: Reduced Chi-squared=%#.8g (Increased)\n",
                iter, lm->new.rchisq);
      }
/*
 * Record the fact that the new iteration did not generate a new best fit.
 */
      was_best = 0;
      break;
    };
  };
/*
 * Reinstate the original current IF state.
 */
  if(set_cif_state(ob, old_if))
    return endfit(ob, mf, quiet, 1);
  return endfit(ob, mf, quiet, 0);
}

/*.......................................................................
 * Private return function of fituvmodel() to perform cleanup.
 */
static int endfit(Observation *ob, Modfit *mf, int quiet, int retcode)
{
/*
 * Record the estimated errors on the fitted parameters.
 */
  if(mf && update_model_errors(mf))
    retcode = 1;
/*
 * Show the fitted model.
 */
  if(!quiet && mf && show_model(ob, ob->newmod, stdout))
    retcode = 1;
/*
 * Return any resources used during the fit.
 */
  del_Modfit(mf);
  return retcode;
}

/*.......................................................................
 * Count the number of free parameters in a given model.
 *
 * Input:
 *  mod    Model *   The model to count from.
 * Output:
 *  return   int     The number of free parameters counted.
 */
static int mod_nfree(Model *mod)
{
  int nfree = 0;   /* The number of free parameters */
/*
 * Are there any components?
 */
  if(mod!=NULL && mod->ncmp>0) {
    Modcmp *cmp;
    for(cmp=mod->head; cmp; cmp = cmp->next) {
      int freepar = cmp->freepar;
      if(freepar) {
/*
 * Ensure that delta components can't have any of the major,ratio or phi
 * parameters variable.
 */
	if(cmp->type==M_DELT) {
	  freepar = (cmp->freepar &= ~(M_MAJOR | M_RATIO | M_PHI));
	};
/*
 * Count free parameters.
 */
	nfree += (freepar & M_FLUX ? 1:0) +
	         (freepar & M_CENT ? 2:0);  /* Note M_CENT implies X and Y */
/*
 * Elliptical parameters?
 */
	if(freepar & (M_MAJOR | M_RATIO | M_PHI)) {
/*
 * M_RATIO must be set if X and Y are to be varied. The only case where
 * X and Y are allowed to be fixed is when a circular aspect has been
 * specified, in which case X and Y are both 0.0.
 */
	  if(cmp->ratio != 1.0)
	    cmp->freepar |= M_RATIO;
	  nfree += cmp->freepar & M_RATIO ? 3:1;  /* X,Y,Z or just Z */
	};
/*
 * Spectral index?
 */
	if(cmp->freepar & M_SPCIND)
	  nfree++;
      };
    };
  };
  return nfree;
}

/*.......................................................................
 * Return the current values of the model free parameters in a provided
 * array.
 *
 * Input:
 *  obj    void *  The Modfit descriptor cast to (void *).
 *  nfree   int    The number of free parameters.
 * Input/Output:
 *  pars double *  An array of nfree elements will be provided. On return
 *                 this should contain the values of each free parameter
 *                 in the model.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static GETFREE(getfree)
{
  Modfit *mf = (Modfit *) obj;  /* The model fit descriptor */
  Modcmp *cmp;                  /* The model component being processed */
/*
 * Sanity check nfree.
 */
  if(nfree != mf->nfree) {
    lprintf(stderr, "getfree: Inconsistent number of free parameters.\n");
    return 1;
  };
/*
 * Find the free parameters in mf->mod and assign their values sequentially
 * to pars[].
 */
  for(cmp=mf->mod->head; cmp; cmp = cmp->next) {
    int freepar = cmp->freepar;
    if(freepar) {
      if(freepar & M_FLUX)
	*pars++ = cmp->flux;
      if(freepar & M_CENT) {
	*pars++ = cmp->x * mf->uvrmax;  /* Re-normalized units */
	*pars++ = cmp->y * mf->uvrmax;
      };
/*
 * The elliptical free parameters.
 */
      if(freepar & (M_MAJOR | M_RATIO | M_PHI)) {
	double anorm = cmp->major * mf->uvrmax;
	double half_aa = 0.5 * anorm * anorm;
	double gg = cmp->ratio * cmp->ratio;
/*
 * X and Y are free parameters unless a circular aspect has been requested.
 */
	if(freepar & M_RATIO) {
	  *pars++ = half_aa * (1.0 - gg) * cos(2.0*cmp->phi);  /* X */
	  *pars++ = half_aa * (1.0 - gg) * sin(2.0*cmp->phi);  /* Y */
	};
	*pars++ = half_aa * (1.0 + gg); /* Z */
      };
/*
 * Spectral index.
 */
      if(freepar & M_SPCIND)
	*pars++ = cmp->spcind;
    };
  };
  return 0;
}

/*.......................................................................
 * Set the current values of the model free parameters from a provided
 * array.
 *
 * Input:
 *  obj    void *  The Modfit descriptor cast to (void *).
 *  nfree   int    The number of free parameters.
 *  pars double *  An array of nfree elements will be provided, containing
 *                 the values to assign to each free parameter in the model.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static SETFREE(setfree)
{
  Modfit *mf = (Modfit *) obj;  /* The model fit descriptor */
  Modcmp *cmp;                  /* The model component being processed */
/*
 * Sanity check nfree.
 */
  if(nfree != mf->nfree) {
    lprintf(stderr, "setfree: Inconsistent number of free parameters.\n");
    return 1;
  };
/*
 * Find the free parameters in mf->mod and assign their values sequentially
 * from pars[].
 */
  for(cmp=mf->mod->head; cmp; cmp = cmp->next) {
    int freepar = cmp->freepar;
    if(freepar) {
/*
 * Instate the flux if variable.
 */
      if(freepar & M_FLUX) {
	float flux = *pars++;
	cmp->flux = flux;
      };
/*
 * Instate the new component centroid position if variable.
 */
      if(freepar & M_CENT) {
	cmp->x = *pars++ / mf->uvrmax;
	cmp->y = *pars++ / mf->uvrmax;
      };
/*
 * Instate the elliptical parameters if variable.
 */
      if(freepar & (M_MAJOR | M_RATIO | M_PHI)) {
	double renorm = mf->uvrmax * mf->uvrmax;
/*
 * Are all of the elliptical parameters variable?
 */
	if(freepar & M_RATIO) {
	  double x = *pars++ / renorm;
	  double y = *pars++ / renorm;
	  double z = *pars++ / renorm;
	  double xyrad = sqrt(x*x+y*y);
/*
 * The minor axis extent is given by sqrt(Z - sqrt(X*X+Y*Y)). Thus if
 * Z < sqrt(X*X+Y*Y) both the minor axis and the axial ratio are imaginary.
 * This is clearly not physical. In such cases we will modify Z such
 * that Z==sqrt(X*X+Y*Y). This has no effect on the major axis length or
 * position angle, but has the effect of setting the axial ratio to zero.
 */
	  if(z < xyrad)
	    z = xyrad;
/*
 * Compute the output model-component parameters.
 */
	  cmp->major = sqrt(fabs(z + xyrad));
	  if(cmp->major == 0.0)
	    cmp->ratio = 1.0;
	  else
	    cmp->ratio = sqrt(fabs(z - xyrad)) / cmp->major;
	  cmp->phi = x==0.0 && y==0.0 ? 0.0 : 0.5 * atan2(y,x);
/*
 * Handle the circular case where only Z is varied.
 */
	} else {
	  double z = *pars++ / renorm;
	  cmp->major = sqrt(fabs(z));
	};
      };
/*
 * Spectral index.
 */
      if(freepar & M_SPCIND)
	cmp->spcind = *pars++;
    };
  };
  return 0;
}

/*.......................................................................
 * Return pertinent details of the next visibility datum to be processed.
 *
 * Input:
 *  obj      void *   The Modfit descriptor.
 *  nfree     int     The number of free parameters in the model.
 * Input/Output:
 *  dy     double *   On output *dy will be assigned with the model vs.
 *                    observed real or imaginary part of the visibility
 *                    residual as (data - model).
 *  wt     double *   On output *wt will be assigned with the weight of
 *                    the data point (1/sigma^2).
 *  mgrad  double *   An array of nfree elements. On output this will be
 *                    assigned with the real or imaginary part of the
 *                    partial derivative of the model wrt each of the
 *                    free parameters.
 * Output:
 *  return    int     0 - End of data.
 *                    1 - Valid datum returned.
 *                   -1 - Error.
 */
static GETNEXT(getnext)
{
  Modfit *mf = (Modfit *) obj; /* The modelfit resource container */
  int i;
/*
 * Is the latest visibility complete?
 */
  if(mf->done) {
    Observation *ob = mf->ob;
    Visibility *vis;
    float uu,vv;       /* The UV coordinates of the visibility */
    float uvrad;       /* The UV radius of the visibility */
    int skip;          /* True if the current visibility isn't usable */
/*
 * Did the previous call to this function process the last visibility
 * datum?
 */
    if(mf->eod) {
      mf->eod = 0;  /* Enable the next cycle */
      return 0;     /* Return End-of-data flag */
    };
/*
 * Get the first IF to be processed?
 */
    if(mf->cif < 0) {
      mf->cif = nextIF(ob, 0, 1, 1);
      if(mf->cif < 0 || getIF(ob, mf->cif)) {
	lprintf(stderr,
	 "modfit: Unable to find any IFs that contain selected channels.\n");
	return -1;
      };
    };
/*
 * Find the next unflagged visibility and record the indexes of the
 * visibility that follows it.
 */
    do {
/*
 * Get the latest visibility.
 */
      vis = &ob->sub[mf->isub].integ[mf->itime].vis[mf->ibase];
/*
 * Is the visibility within the specified UV range?
 */
      uu = vis->u * ob->stream.uvscale;
      vv = vis->v * ob->stream.uvscale;
      uvrad = sqrt(uu*uu+vv*vv);
/*
 * Note that in the following loop if we move on to the next IF,
 * then 'vis' will refer to a different visibility, so extract
 * all information from the lastest visibility now.
 *
 * Should we skip this visibility?
 */
      skip = vis->bad || uvrad < mf->uvrmin || uvrad > mf->uvrmax;
/*
 * Compute the UV representation of the model and its derivatives at
 * the UV coordinates of the latest visibility.
 */
      if(!skip && getmodvis(mf, vis))
	return -1;
/*
 * Increment indexes for the next visibility.
 */
      if(++mf->ibase >= ob->sub[mf->isub].nbase) {
	mf->ibase = 0;
	if(++mf->itime >= ob->sub[mf->isub].ntime) {
	  mf->itime = 0;
	  if(++mf->isub >= ob->nsub) {
	    mf->isub = 0;
	    if((mf->cif=nextIF(ob, mf->cif+1, 1, 1)) < 0) {
	      mf->eod = 1; /* End-of-data now, or on the next call. */
	    } else {
	      if(getIF(ob, mf->cif))
		return -1;
	    };
	  };
	};
      };
    } while(skip && !mf->eod);
/*
 * No more usable visibilities to be processed?
 */
    if(skip) {
      mf->eod = 0;        /* Enable the next cycle */
      return 0;           /* Return End-Of-Data code */
    };
  };
/*
 * Return the real part of the latest visibility parameterization?
 */
  if(mf->done) {
    *dy = mf->re;  /* Real (data - model) residual */
    *wt = mf->wt;  /* Weight given to real part of data */
/*
 * Copy the real parts of the model gradient wrt chi-squared into the
 * return array.
 */
    for(i=0; i<mf->nfree; i++)
      mgrad[i] = mf->vp[i].re;
  }
/*
 * Return the imaginary part of the latest visibility parameterization.
 */
  else {
    *dy = mf->im;  /* Imaginary (data - model) residual */
    *wt = mf->wt;  /* Weight given to imaginary part of data */
/*
 * Copy the imaginary parts of the model gradient wrt chi-squared into the
 * return array.
 */
    for(i=0; i<mf->nfree; i++)
      mgrad[i] = mf->vp[i].im;
  };
/*
 * Toggle the visibility status.
 */
  mf->done = !mf->done;
/*
 * Signal a good data point.
 */
  return 1;
}

/*.......................................................................
 * Construct a UV modelfit object.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  mod        Model *  The model containing free parameters to be fitted.
 *  uvmin      float    The minimum UV radius to take visibilities from.
 *  uvmax      float    The maximum UV radius to take visibilities from.
 * Output:
 *  return    Modfit *  The Modfit object, or NULL on error.
 */
static Modfit *new_Modfit(Observation *ob, Model *mod,
			  float uvmin, float uvmax)
{
  Modfit *mf; /* The new model-fit container */
  int nfree;  /* The number of free parameters in 'mod' */
/*
 * Count the number of free parameters in 'mod'.
 */
  if(mod==NULL || (nfree = mod_nfree(mod)) < 1) {
    lprintf(stderr, "fituvmodel: There are no free parameters to be fitted.\n");
    return NULL;
  };
/*
 * Allocate the model-fit container.
 */
  mf = (Modfit *) malloc(sizeof(Modfit));
  if(mf==NULL) {
    lprintf(stderr, "Insufficient memory to model fit.\n");
    return NULL;
  };
/*
 * Initialize it at least to the point at which it can safely be handed to
 * del_Modfit().
 */
  mf->lm = NULL;
  mf->ob = ob;
  mf->mod = mod;
  mf->nfree = nfree;
  mf->cif = -1;   /* This tells getnext() to read the first IF */
  mf->isub = 0;
  mf->itime = 0;
  mf->ibase = 0;
  mf->done = 1;
  mf->eod = 0;
  mf->vp = NULL;
  mf->re = 0.0;
  mf->im = 0.0;
  mf->wt = 0.0;
  mf->uvrmin = 1.0;
  mf->uvrmax = 1.0;
/*
 * Allocate the array for recording complex partial derivatives.
 */
  mf->vp = (Vispar *) malloc(sizeof(Vispar) * nfree);
  if(mf->vp==NULL) {
    lprintf(stderr, "Insufficient memory to model fit.\n");
    return del_Modfit(mf);
  };
/*
 * Determine the maximum UV radius in the observation. This will
 * be used to scale down U and V, and scale up the major and minor
 * axes, both the avoid overflows in the hessian matrix (where U and
 * V enter to the power 4), and to reduce the range of magnitudes in
 * the Hessian.
 */
  {
    UVrange *uvr = uvrange(ob, 1, 0, uvmin, uvmax);
    if(uvr==NULL)
      return del_Modfit(mf);
    mf->uvrmin = uvr->uvrmin;
    mf->uvrmax = uvr->uvrmax;
  };
/*
 * Get a Levensberg-Marquardt non-linear least-squares object.
 * For safety this has to be done after the rest of 'mf' has been
 * initialized since new_Lmfit() makes no guarantees that the
 * getfree,setfree and getnext functions are not called upon to operate
 * on 'mf' - which must thus be a complete object.
 */
  mf->lm = new_Lmfit((void *)mf, nfree, getfree, setfree, getnext);
  if(mf->lm==NULL)
    return del_Modfit(mf);
/*
 * Return the new descriptor.
 */
  return mf;
}

/*.......................................................................
 * Delete a Modfit descriptor.
 *
 * Input:
 *  mf     Modfit *   The descriptor to be deleted.
 * Output:
 *  return Modfit *   Allways NULL. Use like, mf=del_Modfit(mf);
 */
static Modfit *del_Modfit(Modfit *mf)
{
  if(mf) {
    if(mf->vp)
      free(mf->vp);
    if(mf->lm)
      del_Lmfit(mf->lm);
  };
  return NULL;
}

/*.......................................................................
 * Calculate the UV plane complex representation of the current variable
 * model and its derivatives wrt each of the free parameters.
 *
 * Input:
 *  mf      Modfit *  The model-fit descriptor.
 *  vis Visibility *  The descriptor of the visibility.
 * Output:
 *  mf->vp[0..nfree-1] The complex partial derivatives of the model wrt
 *                     each of the model free parameters.
 *  mf->re             The real part of the (data - model) residual
 *                     where the model includes both the variable model
 *                     mf->mod and the established model.
 *  mf->im             The imaginary part of the (data - model) residual
 *                     where the model includes both the variable model
 *                     mf->mod and the established model.
 *  mf->wt             The weight of the real and imaginary parts of the
 *                     visibility.
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int getmodvis(Modfit *mf, Visibility *vis)
{
  Vispar *vp = mf->vp;   /* The model vs. free-parameter partial derivatives */
  UVstream *uvs = &mf->ob->stream;   /* Stream-specific parameters */
  double uu = vis->u * uvs->uvscale; /* Visibility U in wavelengths */
  double vv = vis->v * uvs->uvscale; /* Visibility V in wavelengths */
  double uun = uu / mf->uvrmax;      /* Re-normalized version of uu */
  double vvn = vv / mf->uvrmax;      /* Re-normalized version of vv */
  double freq = getfreq(mf->ob, mf->cif); /* The frequency of the current IF */
  Modcmp *cmp;                       /* The model component being processed */
  int i;
/*
 * Get the current (data - established_model) data residual and weight.
 */
  mf->re = vis->amp * cos(vis->phs) - vis->modamp * cos(vis->modphs);
  mf->im = vis->amp * sin(vis->phs) - vis->modamp * sin(vis->modphs);
  mf->wt = vis->wt;
/*
 * Clear the partial-derivative output array.
 */
  for(i=0; i<mf->nfree; i++)
    mf->vp[i].re = mf->vp[i].im = 0.0;
/*
 * Loop through the components of the variable model.
 */
  for(cmp=mf->mod->head; cmp; cmp=cmp->next) {
/*
 * Since all model component types are even functions, the only
 * contribution to the model visibility phase is from the centroid
 * position of the component.
 */
    double cmpphs = twopi * (uu*cmp->x + vv*cmp->y); /* Component phase */
    double cmpamp;                                   /* Component amplitude */
    double cmpre;   /* Real part of model-component visibility */
    double cmpim;   /* Imaginary part of model-component visibility */
/*
 * Pre-compute useful constants.
 */
    double sinphi = sin(cmp->phi);
    double cosphi = cos(cmp->phi);
/*
 * Compute the elliptically stretched UV radius (also scaled by pi * major
 * axis for convenience).
 */
    double tmpa = (uu * cosphi - vv * sinphi) * cmp->ratio;
    double tmpb = (uu * sinphi + vv * cosphi);
    double uvrad = pi * cmp->major * sqrt(tmpa*tmpa + tmpb*tmpb);
/*
 * Get the spectral-index scale factor.
 */
    double si = cmp->spcind==0.0 ? 1.0 : pow(freq/cmp->freq0, cmp->spcind);
/*
 * Get the primary beam scale factor.
 */
    double pb = pb_bl_factor(mf->ob->sub + mf->isub, mf->ibase, freq,
			     calc_pointing_offset(mf->ob, cmp->x, cmp->y));
/*
 * Get the potentially frequency dependent flux of the component.
 */
    double flux = cmp->flux * si * pb;
/*
 * Get the bitmap of free-parameter designations.
 */
    int freepar = cmp->freepar;
/*
 * Limit uvrad to sensible values to prevent underflow,overflow and
 * divide-by-zero errors.
 */
    if(uvrad < 1.0e-9)
      uvrad = 1.0e-9;
/*
 * Compute model-component visibility amplitude.
 */
    switch(cmp->type) {
    case M_DELT:
      cmpamp = flux;
      break;
    case M_GAUS:
      cmpamp = flux * (uvrad < 12.0 ? exp(-0.3606737602 * uvrad*uvrad) : 0.0);
      break;
    case M_DISK:
      cmpamp = 2.0f * flux * c_besj1(uvrad)/uvrad;
      break;
    case M_ELLI:
      cmpamp = 3.0f * flux * (sin(uvrad)-uvrad*cos(uvrad))/(uvrad*uvrad*uvrad);
      break;
    case M_RING:
      cmpamp = flux * c_besj0(uvrad);
      break;
    case M_RECT:
      lprintf(stderr, "getmodvis: Rectangular components are not supported.\n");
      return 1;
      break;
    case M_SZ:
      cmpamp = flux * (uvrad < 50.0 ? exp(-uvrad) : 0.0) / uvrad;
      break;
    default:
      lprintf(stderr, "Unknown model component type: %d\n", cmp->type);
      return 1;
      break;
    };
/*
 * Compute the real and imaginary parts of the model-component
 * visibility.
 */
    cmpre = cmpamp * cos(cmpphs);
    cmpim = cmpamp * sin(cmpphs);
/*
 * Now compute partial derivatives for free parameters in the current component.
 *
 * The partial derivative of a component wrt flux.
 */
    if(freepar & M_FLUX) {
      vp->re = cmpre / cmp->flux;
      vp->im = cmpim / cmp->flux;
      vp++;
    };
/*
 * The partial derivatives of a component wrt centroid X and Y positions
 * position.
 */
    if(freepar & M_CENT) {
      vp->re = twopi * uun * -cmpim; /* Derivative wrt X */
      vp->im = twopi * uun * cmpre;
      vp++;
      vp->re = twopi * vvn * -cmpim; /* Derivative wrt Y */
      vp->im = twopi * vvn * cmpre;
      vp++;
    };
/*
 * The elliptical parameterization partial derivatives all require a
 * common factor to be computed, that depends on the component type.
 */
    if(freepar & (M_MAJOR | M_RATIO | M_PHI)) {
      double comfac=0.0;  /* The common factor */
/*
 * Determine the common factor used in the partial derivatives below.
 */
      switch(cmp->type) {
      case M_DELT:
	comfac = 0.0; /* No dependence on major,ratio or phi */
	break;
      case M_GAUS:
	comfac = -0.7213475204 * uvrad;
	break;
      case M_DISK:
	comfac = -2.0 * c_besj2(uvrad)/uvrad;
	break;
      case M_ELLI:
	comfac = (9.0 * cos(uvrad) / uvrad -
		  9.0 * sin(uvrad) / (uvrad * uvrad) +
		  3.0 * sin(uvrad)) / uvrad / uvrad;
	break;
      case M_RING:
	comfac = -c_besj1(uvrad);
	break;
      case M_RECT:
	lprintf(stderr, "modfit: Rectangular components are not supported.\n");
	return 1;
	break;
      case M_SZ:
	comfac = -(uvrad<50.0 ? exp(-uvrad):0.0) * (uvrad+1.0) / uvrad / uvrad;
	break;
      default:
	lprintf(stderr, "modfit: Unknown model component type: %d\n",cmp->type);
	return 1;
	break;
      };
/*
 * Compute useful parameters.
 */
      {
	double newfac = comfac * 0.5 * pi * pi / uvrad;
/*
 * X and Y are free parameters unless a circular aspect has been requested.
 */
	if(freepar & M_RATIO) {
	  double tmp;
/*
 * Partial derivative wrt X.
 */
	  tmp = newfac * (vvn - uun) * (vvn + uun);
	  vp->re = cmpre * tmp;
	  vp->im = cmpim * tmp;
	  vp++;
/*
 * Partial derivative wrt Y.
 */
	  tmp = newfac * (2.0 * uun * vvn);
	  vp->re = cmpre * tmp;
	  vp->im = cmpim * tmp;
	  vp++;
	};
/*
 * Partial derivative wrt Z.
 */
	{
	  double tmp = newfac * (vvn*vvn + uun*uun);
	  vp->re = cmpre * tmp;
	  vp->im = cmpim * tmp;
	  vp++;
	};
      };
    };
/*
 * The partial derivatives of the component wrt the spectral
 * index.
 */
    if(freepar & M_SPCIND) {
      double factor = log(freq/cmp->freq0);
      vp->re = cmpre * factor;
      vp->im = cmpim * factor;
      vp++;
    };
/*
 * Accumulate the residual visibility.
 */
    mf->re -= cmpre;
    mf->im -= cmpim;
  };
/*
 * Job completed succesfully.
 */
  return 0;
}

/*.......................................................................
 * Update the estimated errors of the free parameters of the fit.
 *
 * Input:
 *  mf     Modfit *   The model fitting resource object.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int update_model_errors(Modfit *mf)
{
  Modcmp *cmp;     /* The model component being processed */
  double **covar;  /* The covariance matrix of the fit */
  double *pars;    /* The array of best fit parameter values */
  int npar = 0;    /* The number of parameters processed so far */
  double var;      /* The estimated variance of a paremeter */
  const double tiny = 1e-50;
/*
 * Check the arguments.
 */
  if(!mf) {
    lprintf(stderr, "update_model_errors: Missing Modfit object.\n");
    return 1;
  }
/*
 * Get the covariance matrix of the fit.
 */
  covar = lm_covar(mf->lm);
  if(!covar)
    return 1;
/*
 * Get the best-fit values.
 */
  pars = mf->lm->best.pars;
/*
 * Find the free parameters in mf->mod and assign their errors sequentially
 * from the diagonal of the covariance matrix.
 */
  for(cmp=mf->mod->head; cmp; cmp = cmp->next) {
    int freepar = cmp->freepar;
/*
 * Start with all the uncertainties reset to zero.
 */
    cmp->flux_err = 0.0;
    cmp->x_err = 0.0;
    cmp->y_err = 0.0;
    cmp->major_err = 0.0;
    cmp->ratio_err = 0.0;
    cmp->phi_err = 0.0;
    cmp->spcind_err = 0.0;
/*
 * If any parameters of this component were free to vary, update their
 * statistical uncertainties.
 */
    if(freepar) {
/*
 * Compute the uncertainty in the flux, if variable.
 */
      if(freepar & M_FLUX) {
        var = covar[npar][npar]; /* Get flux variance from covariance diagonal */
        if(var > 0.0) cmp->flux_err = sqrt(var);
        npar++;
      }
/*
 * Compute the uncertainties of component's centroid, if variable.
 */
      if(freepar & M_CENT) {
/*
 * The x-axis position.
 */
        var = covar[npar][npar];  /* Get x variance from covariance diagonal */
        if(var > 0.0) cmp->x_err = sqrt(var) / mf->uvrmax;
        npar++;
/*
 * The y-axis position.
 */
        var = covar[npar][npar];  /* Get y variance from covariance diagonal */
        if(var > 0.0) cmp->y_err = sqrt(var) / mf->uvrmax;
        npar++;
      };
/*
 * Compute errors for the elliptical parameters if variable.
 */
      if(freepar & (M_MAJOR | M_RATIO | M_PHI)) {
/*
 * Get the magnitude renormalization constant.
 */
	double renorm = mf->uvrmax * mf->uvrmax;
        double renorm_sq = renorm * renorm;
/*
 * Start with all of the uncertainties marked as zero. We will
 * override the ones that have known uncertainties below.
 */
        double major_var = 0.0;
        double ratio_var = 0.0;
        double phi_var = 0.0;
/*
 * Are all the elliptical parameters variable?
 */
	if(freepar & M_RATIO) {
/*
 * The intermediate x parameter of the ellipse.
 */
	  double x = pars[npar] / renorm;
          double xx = x*x;
          double x_var = covar[npar][npar] / renorm_sq;
          npar++;
/*
 * The intermediate y parameter of the ellipse.
 */
	  double y = pars[npar] / renorm;
          double yy = y*y;
          double y_var = covar[npar][npar] / renorm_sq;
          npar++;
/*
 * The intermediate z parameter of the ellipse.
 */
	  double z = pars[npar] / renorm;
          double z_var = covar[npar][npar] / renorm_sq;
          npar++;
/*
 * Get the radius of the vector (x,y).
 */
          double rr = xx + yy;
	  double r = sqrt(rr);
/*
 * Temporary variables for holding the fitted major axis dimension.
 */
          double major;
/*
 * Compute the variance of r.
 */
          if(rr > tiny) {
/*
 * Determine the uncertainty in r.
 */
            double r_var = (xx * x_var + yy * y_var) / rr;
/*
 * The minor axis extent is given by sqrt(Z - sqrt(X*X+Y*Y)). Thus if
 * Z < sqrt(X*X+Y*Y) both the minor axis and the axial ratio are
 * imaginary.  This is clearly not physical. In such cases we will
 * modify Z such that Z==sqrt(X*X+Y*Y), which sets the axial ratio
 * (and thus the minor extent) to zero, increases the major axis
 * length, but leaves the position angle unchanged.
 */
            if(z < r) {
              z = r;
              z_var = r_var;
            }
/*
 * Compute the errors on the elliptical dimensions.
 */
            major = sqrt(fabs(z + r));
            if(major > tiny) {
              double minor = sqrt(fabs(z - r));
              double major_sqr = major * major;
              major_var = 0.25 / major_sqr * (z_var + r_var);
              if(minor > tiny) {
                double minor_sqr = minor * minor;
                double minor_var = 0.25 / minor_sqr * (z_var + r_var);
                ratio_var = minor_var / major_sqr
                  + major_var * minor_sqr / major_sqr / major_sqr;
              }
            }
/*
 * Compute the error on the elliptical angle: phi = 0.5 * atan2(y,x)
 */
            if(fabs(x) > tiny) {
              double q = y / x;
              double qq = q*q;
              double tmp = (x * (1.0 + qq));
              phi_var = 0.25 * (y_var + qq * x_var) / tmp / tmp;
            }
          }
/*
 * Handle the circular case where only Z is varied.
 */
	} else {
	  double z = pars[npar] / renorm;
          double z_var = covar[npar][npar] / renorm_sq;
          double major = sqrt(fabs(z));
          npar++;
          if(major > tiny) {
            major_var = 0.25 / major / major * z_var;
          }
	};
/*
 * Record the standard deviations.
 */
        if(major_var > 0.0) cmp->major_err = sqrt(major_var);
        if(ratio_var > 0.0) cmp->ratio_err = sqrt(ratio_var);
        if(phi_var > 0.0) cmp->phi_err = sqrt(phi_var);
      };
/*
 * Spectral index.
 */
      if(freepar & M_SPCIND) {
        double var = covar[npar][npar];
        npar++;
        if(var > 0.0) cmp->spcind_err = sqrt(var);
      }
    };
  };
  return 0;
}

/*.......................................................................
 * Write the details of a fitted model to a file or the terminal.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  mod        Model *  The model containing free parameters to be fitted.
 *  fd          FILE *  The file descriptor of an open file. I am assuming
 *                      that all users of this function will want to start
 *                      the file with a descriptive header before calling
 *                      this routine - hence the requirement for an open
 *                      file rather than the name of a file to be opened.
 * Output:
 *   return   int   0 - No write errors.
 */
int show_model(Observation *ob, Model *mod, FILE *fd)
{
  Modcmp *cmp;   /* The component being described */
/*
 * Check the arguments.
 */
  if(!ob || !mod) {
    lprintf(stderr, "show_model: Missing argument(s).\n");
    return 1;
  }
/*
 * If no output file was specified, substitude stdout.
 */
  if(!fd)
    fd = stdout;
/*
 * Label the file columns.
 */
  lprintf(fd, "#    Flux (Jy)          East (arcsec)        North (arcsec)    Shape   R.A. (deg)       Dec (deg)  Major FWHM (arcsec) Minor FWHM (arcsec)  Theta (deg)     Freq (Hz)      Spectral Index  \n");
  lprintf(fd, "#  Value     Stdev      Value     Stdev      Value     Stdev                                          Value    Stdev     Value    Stdev     Value  Stdev                   Value     Stdev \n");
  lprintf(fd, "#------------------  -------------------  ------------------- ------  -------------- -------------- -----------------  -----------------  ---------------  -----------  -------------------\n");

/*
 * No components to be written?
 */
  if(mod==NULL || mod->ncmp==0)
    return 0;
/*
 * Write the components according to their types.
 */
  for(cmp=mod->head; cmp != NULL;cmp=cmp->next) {
/*
 * Compute the RA, Dec of the component’s centroid.
 */
    double ra = lmtora(ob->source.ra, ob->source.dec,
                       -ob->geom.east + cmp->x,
                       -ob->geom.north + cmp->y, ob->proj);
    double dec = lmtodec(ob->source.ra, ob->source.dec,
                         -ob->geom.east + cmp->x,
                         -ob->geom.north + cmp->y, ob->proj);
/*
 * All component have a flux and an x,y position.
 */
    lprintf(fd,"% 11.4e %7.1e", cmp->flux, cmp->flux_err);
    lprintf(fd,"  % 11.4e %7.1e", cmp->x * rtoas, cmp->x_err * rtoas);
    lprintf(fd,"  % 11.4e %7.1e", cmp->y * rtoas, cmp->y_err * rtoas);
/*
 * Show the component type.
 */
    lprintf(fd," %6s", modtyp_name(cmp->type));
/*
 * Show the Right Ascension and Declination of the above centroid, giving
 * both fields a precision of one micro-arcsecond.
 */
    lprintf(fd,"  %14.10f % 14.10f", ra * rtod, dec * rtod);
/*
 * All component types except delta functions have an elliptical
 * shape.
 */
    if(cmp->type != M_DELT || cmp->freq0 > 0.0) {
      double minor = cmp->major * cmp->ratio;
      double minor_err = cmp->major * cmp->ratio_err;
      lprintf(fd," %9.3e %7.1e", cmp->major * rtoas,
              cmp->major_err * rtoas);
      lprintf(fd,"  %9.3e %7.1e", minor * rtoas,
              minor_err * rtoas);
      lprintf(fd,"  %# 7.2f %7.1e", cmp->phi * rtod,
              cmp->phi_err * rtod);
      if(cmp->freq0 > 0.0) {
        lprintf(fd, "  %#11.6g", cmp->freq0);
        lprintf(fd, "  % 10.3e  %7.1e", cmp->spcind, cmp->spcind_err);
      };
    };
/*
 * New line.
 */
    lputc('\n', fd);
/*
 * Check for write errors.
 */
    if(ferror(fd)) {
      lprintf(stderr, "Error writing model file\n");
      return -1;
    };
  };
/*
 * Finished successfully.
 */
  return 0;
}

