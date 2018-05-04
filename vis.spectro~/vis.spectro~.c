/*
	vis.spectro~ -- new and improved spectroscope
 
    Part of a new group of smart audio visualization objects for Max which includes an oscilloscope,
    a spectroscope, a 3D spectroscope, smart meters, and maybe more
*/

#include "ext.h"
#include "ext_obex.h"
#include "jpatcher_api.h"
#include "jgraphics.h"
#include "z_dsp.h"
#include "fft_mayer.proto.h"


#define FFT_DEFAULT_POINTS 2048
#define FFT_MAX_POINTS	4096
#define FFT_MIN_POINTS	16


typedef struct _spectro
{
    
	t_pxjbox u_box;
	void *u_out;
	t_jrgba u_outline;
	t_jrgba u_background;
    t_jrgba u_samples;
    
    double u_bordersize;
    double u_linewidth;
    int u_gridwidth;
    double u_gridheight;
    
    long		f_inverse;
    long		f_fftsize;		// size
    long		f_interval;		// sample count before doing another FFT (i.e. hop size, but always >= to fftsize)
    long		f_phase;		// phase of FFT start (relative to global counter)
    long		f_countdown;	// vector count for phase offset fft
    double		f_1overpts;		// (1.0/fftsize) for inverse FFT scaling
    
    t_double	*f_realin;		// where transform is done
    t_double	*f_imagin;		// where transform is done
    t_double	*f_realout;		// where transform is done
    t_double	*f_imagout;		// where transform is done
    
    t_double	*f_realinptr;	// input fill ptr into realin
    t_double	*f_imaginptr;	// input fill ptr into imagin
    t_double	*f_realoutptr;	// output ptr into realout
    t_double	*f_imagoutptr;	// output ptr into imagout
    
} t_spectro;


void *spectro_new(t_symbol *s, long argc, t_atom *argv);
void spectro_free(t_spectro *x);
void spectro_assist(t_spectro *x, void *b, long m, long a, char *s);
void spectro_paint(t_spectro *x, t_object *patcherview);
void spectro_getdrawparams(t_spectro *x, t_object *patcherview, t_jboxdrawparams *params);

void spectro_setphase(t_spectro *x, long phase);
void spectro_realimag_perform64(t_spectro *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void spectro_dsp64(t_spectro *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);


static t_class *s_spectro_class;

static t_symbol	*ps_fft;
static t_symbol	*ps_ifft;

void ext_main(void *r)
{
	t_class *c;

	c = class_new("vis.spectro~", (method)spectro_new, (method)spectro_free, sizeof(t_spectro), 0L, A_GIMME, 0);

	c->c_flags |= CLASS_FLAG_NEWDICTIONARY;
	jbox_initclass(c, JBOX_FIXWIDTH | JBOX_COLOR);
    class_dspinitjbox(c);

	class_addmethod(c, (method)spectro_paint,		"paint",	A_CANT, 0);
	class_addmethod(c, (method)spectro_assist,		"assist",	A_CANT, 0);
    class_addmethod(c, (method)spectro_dsp64,		"dsp64", A_CANT, 0);

    
    
	// attributes
	CLASS_STICKY_ATTR(c, "category", 0, "Color");
    
	CLASS_ATTR_RGBA(c, "bgcolor", 0, t_spectro, u_background);
    CLASS_ATTR_BASIC(c, "bgcolor", 0);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bgcolor", 0, "1. 1. 1. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bgcolor",0,"rgba","Background Color");
    
    CLASS_ATTR_RGBA(c, "samplecolor", 0, t_spectro, u_samples);
    CLASS_ATTR_BASIC(c, "samplecolor", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "samplecolor", 0, "0. 0. 255. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"samplecolor",0,"rgba","Sample Color");

	CLASS_ATTR_RGBA(c, "bordercolor", 0, t_spectro, u_outline);
    CLASS_ATTR_BASIC(c, "bordercolor", 0);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bordercolor", 0, "0. 0. 0. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bordercolor",0,"rgba","Border Color");

	CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    CLASS_STICKY_ATTR(c, "category", 0, "Dimensions");
    
    CLASS_ATTR_DOUBLE(c, "bordersize", 0, t_spectro, u_bordersize);
    CLASS_ATTR_BASIC(c, "bordersize", 0);
    CLASS_ATTR_DEFAULTNAME(c, "bordersize", 0, "1.");
    CLASS_ATTR_STYLE_LABEL(c, "bordersize", 0, "double", "Border Size");
    
    CLASS_ATTR_DOUBLE(c, "linewidth", 0, t_spectro, u_linewidth);
    CLASS_ATTR_BASIC(c, "linewidth", 0);
    CLASS_ATTR_DEFAULTNAME(c, "linewidth", 0, "3.");
    CLASS_ATTR_STYLE_LABEL(c, "linewidth", 0, "double", "Line Width");
    
    CLASS_STICKY_ATTR_CLEAR(c, "category");

	CLASS_ATTR_DEFAULT(c,"patching_rect",0, "0. 0. 450. 200.");

	class_register(CLASS_BOX, c);
	s_spectro_class = c;
}

void spectro_assist(t_spectro *x, void *b, long m, long a, char *s)
{
	if (m == 1)		//inlet
		sprintf(s, "(signal) Audio Input");
}

void spectro_paint(t_spectro *x, t_object *patcherview)
{
	t_rect rect;
	t_jgraphics *g = (t_jgraphics *) patcherview_get_jgraphics(patcherview);		// obtain graphics context
    t_jgraphics *h = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
	jbox_get_rect_for_view((t_object *)x, patcherview, &rect);


	// paint border
	jgraphics_set_source_jrgba(g, &x->u_outline);
	jgraphics_set_line_width(g, x->u_bordersize);
	jgraphics_rectangle(g, x->u_bordersize / 2,
                        x->u_bordersize / 2,
                        rect.width - x->u_bordersize,
                        rect.height - x->u_bordersize);
    
	jgraphics_stroke(g);
    
    //paint background
    jgraphics_set_source_jrgba(g, &x->u_background);
    jgraphics_rectangle(g, x->u_bordersize / 2,
                        x->u_bordersize / 2,
                        rect.width - x->u_bordersize,
                        rect.height - x->u_bordersize);
    
    jgraphics_fill(g);
    
    
	// paint samples
    jgraphics_set_line_width(h, x->u_linewidth);
    jgraphics_set_source_jrgba(h, &x->u_samples);
    
    
    //prep dimensions for waveform drawing
    x->u_gridwidth = rect.width - x->u_bordersize;
    x->u_gridheight = rect.height - x->u_bordersize;

    
    //draw waveform of zero if audio is off
    if(sys_getdspstate() == 0){
        for(int i = 0; i < x->u_gridwidth; i++){
            if(i == 0){
                jgraphics_move_to(g, x->u_bordersize / 2, x->u_gridheight);
                i++;
            } else {
                jgraphics_line_to(g, (double) i + (x->u_bordersize / 2), x->u_gridheight);
                i++;
            }
        }
        
    } else {
        
        //rms scaling
        /*float rms = 0;
        for (int i = 0; i < x->f_fftsize; i++){
            //magnitute calculation
            rms += x->f_realout[i] * x->f_imagout[i];
        }
        rms = sqrt(rms / x->f_fftsize);
        float scale = 1 / rms;
        for (int i = 0; i < x->f_fftsize; i++){
            x->f_realout[i] *= scale;
            x->f_imagout[i] *= scale;
        }
        */
        
        double sample[x->f_fftsize];

        //draw spectroscope
        for(int i = 0; i < x->f_fftsize / 2; i++){
            //scale samples to current gridheight
            float px = ((float) i / (x->f_fftsize / 2)) * x->u_gridwidth;
            
            sample[i] = (sqrt(pow(x->f_imagout[i], 2)+pow(x->f_realout[i], 2)));
            if(i == 0){
                jgraphics_move_to(h, x->u_bordersize / 2, x->u_gridheight - (x->u_bordersize / 2) - sample[i]);
            } else {
                jgraphics_line_to(h, (float) px + (x->u_bordersize / 2), x->u_gridheight - sample[i]);
            }
        }
    }
    
    jgraphics_set_line_width(h, x->u_linewidth);
    jgraphics_stroke(h);
}


void spectro_getdrawparams(t_spectro *x, t_object *patcherview, t_jboxdrawparams *params)
{
	params->d_bordercolor.alpha = 0;
	params->d_boxfillcolor.alpha = 0;
}

void spectro_free(t_spectro *x)
{
    dsp_freejbox((t_pxjbox *)x);
	jbox_free((t_jbox *)x);
}

void *spectro_new(t_symbol *s, long argc, t_atom *argv)
{
	t_spectro *x = NULL;
	t_dictionary *d = NULL;
	long boxflags;
    
    long fftsize = FFT_DEFAULT_POINTS;
    long interval = fftsize;
    long phase = 0;
    long inverse = (s == ps_ifft);
    
    if (!(d = object_dictionaryarg(argc,argv)))
        return NULL;
    
    x = (t_spectro *)object_alloc(s_spectro_class);
    
    x->f_fftsize = FFT_DEFAULT_POINTS; // FFT size
    x->f_interval = interval;
    x->f_phase = phase; // typed-in 3rd arg
    x->f_inverse = inverse;
    
    x->f_realin = 0L; // init pointers just to be safe
    x->f_imagin = 0L;
    x->f_realout = 0L;
    x->f_imagout = 0L;
    x->f_realin = (t_double *)sysmem_newptr(fftsize * 4 * sizeof(t_double));
    x->f_imagin = x->f_realin + fftsize;
    x->f_realout = x->f_realin + (fftsize*2);
    x->f_imagout = x->f_realin + (fftsize*3);
    
    x->f_realinptr = x->f_realin;
    x->f_realoutptr = x->f_realout;
    x->f_imaginptr = x->f_imagin;
    x->f_imagoutptr = x->f_imagout;
    x->f_countdown = 0;
    x->f_1overpts = 1. / x->f_fftsize;
    
    
	boxflags = 0
			   | JBOX_DRAWFIRSTIN
			   | JBOX_NODRAWBOX
			   | JBOX_DRAWINLAST
			   | JBOX_TRANSPARENT
               | JBOX_GROWBOTH
			   | JBOX_DRAWBACKGROUND
			   ;
    
	jbox_new((t_jbox *)x, boxflags, argc, argv);
    x->u_box.z_box.b_firstin = (void *)x;
    dsp_setupjbox((t_pxjbox *)x,2);
    outlet_new(x, "signal");
    outlet_new(x, "signal");
    outlet_new(x, "signal");
	attr_dictionary_process(x,d);
	jbox_ready((t_jbox *)x);
	return x;
}

void spectro_setphase(t_spectro *x, long phase)
{
    x->f_phase = phase;
}


void spectro_dsp64(t_spectro *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    long phase; // typed-in 3rd arg
    long blocksize = maxvectorsize;
    
    x->f_countdown = 0;
    x->f_realinptr = x->f_realin;
    x->f_realoutptr = x->f_realout;
    x->f_imaginptr = x->f_imagin;
    x->f_imagoutptr = x->f_imagout;
    
    // adjust phase based on interval and blocksize
    phase = x->f_phase % x->f_interval;
    if (phase % blocksize) {
        if (phase > blocksize)
            phase -= (phase % blocksize);
        else
            phase = blocksize;
        phase %= x->f_interval; // just in case blocksize >= interval
        object_error((t_object *)x, "%s: phase must be multiple of %ld, setting to %ld",x->f_inverse?ps_ifft->s_name:ps_fft->s_name, blocksize,phase);
    }
    
    if (phase < x->f_fftsize) {
        x->f_realoutptr = x->f_realout + phase;
        x->f_imagoutptr = x->f_imagout + phase;
    }
    else {
        x->f_realoutptr = x->f_realout + x->f_fftsize;
        x->f_imagoutptr = x->f_imagout + x->f_fftsize;
    }
    
    if (x->f_fftsize + phase >= x->f_interval) {
        x->f_countdown = 0;
        x->f_realinptr = x->f_realin + (x->f_fftsize + phase - x->f_interval);
        x->f_imaginptr = x->f_imagin + (x->f_fftsize + phase - x->f_interval);
    }
    else {
        x->f_countdown = (x->f_interval - phase - x->f_fftsize) / blocksize;
        x->f_realinptr = x->f_realin;
        x->f_imaginptr = x->f_imagin;
    }
    
    // zero buffers in case they have junk in them (i.e. now makes no glitches when dsp is turned on)
    memset(x->f_realin, 0, sizeof(*x->f_realin) * x->f_fftsize*4);
    
    dsp_add64(dsp64, (t_object *)x, (t_perfroutine64)spectro_realimag_perform64, 0, NULL);

}

void spectro_realimag_perform64(t_spectro *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) // based on fft_complex_perform()
{
    t_double	*inreal = ins[0];
    t_double	*inimag = ins[1];
    t_double	*outreal = outs[0];
    t_double	*outimag = outs[1];
    t_double	*outsync = outs[2];
    
    int			n = (int) sampleframes;
    int			m, np, count;
    long		fftsize = x->f_fftsize;
    long		interval = x->f_interval;
    int			inverse = (int) x->f_inverse;
    t_double	*b;
    t_double	*ei = x->f_realin + fftsize;						// just do the real buf as the test
    
    if (n > fftsize) {
        np = (int) fftsize;
        count = n / fftsize;
    }
    else {
        np = n;
        count = 1;
    }
    
    while (count--) {
        m = np;														// we need this for sync outlet
        
        if (x->f_countdown) {										// if phase offset then decrement vector countdown
            x->f_countdown--;
            inreal += np;											// need to do this in case sigvs>fftsize
            inimag += np;
            memset(outreal, 0, sizeof(*outreal) * np);
            memset(outimag, 0, sizeof(*outimag) * np);
            outreal += np;											// need to do this in case sigvs>fftsize
            outimag += np;
            while (m--)
                *outsync++ = 0.0;									// output zero while no fft
        }
        else {														// otherwise move input samples into and out of fft buffers
            double q;
            
            // buffer input to fftinput buffer
            memcpy(x->f_realinptr, inreal, sizeof(*inreal) * np);
            memcpy(x->f_imaginptr, inimag, sizeof(*inreal) * np);
            
            // increment pointers
            x->f_realinptr += np;
            x->f_imaginptr += np;
            inreal += np;											// need to do this in case sigvs>fftsize
            inimag += np;
            b = x->f_realinptr;										// to compare later on and see if we need to do a fft
            
            // buffer fftoutput buffer to output
            memcpy(outreal, x->f_realoutptr, sizeof(*inreal) * np);
            memcpy(outimag, x->f_imagoutptr, sizeof(*inreal) * np);
            
            q = (int)(x->f_realoutptr - x->f_realout);				// where to start sync count
            while (m--) {
                *outsync++ = q;
                q += 1.0;
            }
            
            // increment pointers
            x->f_realoutptr += np;
            x->f_imagoutptr += np;
            outreal += np;											// need to do this in case sigvs>fftsize
            outimag += np;
            
            if (b == ei) {											// if realinptr = realin+fftsize
                // this step could be eliminated if signals were unique
                memcpy(x->f_realout, x->f_realin, sizeof(*x->f_realin) * fftsize*2);
                
                if (inverse) {
                    t_double	*real = x->f_realout;
                    t_double	*imag = x->f_imagout;
                    t_double	mult = x->f_1overpts;
                    int			temp = (int) fftsize;						// so fftsize will always be valid
                    
                    ifft64((int) fftsize, real, imag);
                    while (temp--) {
                        *real++ *= mult;
                        *imag++ *= mult;
                    }
                }
                else {
                    // this part does not follow fft_complex_perform()
                    int			nn;
                    t_double	*rout, *iout;
                    t_double	*real = x->f_realout;
                    t_double	*imag = x->f_imagout;
                    
                    realfft64((int) fftsize, x->f_realout);
                    
                    // re-arrange the real fft into complex fft frames
                    rout = real+fftsize-1;
                    iout = imag+fftsize-1;
                    nn = ((int) fftsize/2)-1;
                    
                    *imag++ = 0.0;									// imag at 0 Hz = 0.
                    //*real++;
                    
                    while (nn--) {
                        *imag++ = -(*rout);
                        *iout-- = *rout;
                        *rout-- = *real++;
                    }
                    *imag++ = 0.0;									// imag at sr/2 Hz = 0.
                }
                
                // reset pointers
                x->f_realoutptr = x->f_realout;
                x->f_imagoutptr = x->f_imagout;
                x->f_realinptr = x->f_realin;
                x->f_imaginptr = x->f_imagin;
                x->f_countdown = (interval - fftsize) / np; // not sure this has always worked correctly - to be verified
            }
        }
    }
    
    jbox_redraw((t_jbox *)x);

}

