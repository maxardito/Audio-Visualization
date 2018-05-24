/*
	vis.spectro~ -- new and improved spectroscope
 
    Max Ardito (2018)
 
    Part of a new group of smart audio visualization objects for Max which includes an oscilloscope,
    a spectroscope, a 3D spectroscope, smart meters, and maybe more
*/

#include "ext.h"
#include "ext_obex.h"
#include "jpatcher_api.h"
#include "jgraphics.h"
#include "z_dsp.h"
#include "fft_mayer.proto.h"

#define FFT_MIN_BINS 128
#define FFT_MAX_BINS 2048

typedef struct bin{
    unsigned int number;
    int freq;
    double amp;
} t_bin;

typedef struct _spectro
{
    
	t_pxjbox u_box;
	void *u_out;
	t_jrgba u_outline;
	t_jrgba u_background;
    t_jrgba u_samples;
    t_jrgba u_text;
    t_jrgba u_verticalLine;
    t_jrgba u_binGrid;
    int u_fontSize;
    long u_ampOn;
    
    t_jtextlayout	*mytxt;
    t_jfont	*myfont;
    t_jrgba textcolor;
    
    double u_mouseX;
    bool u_mouseover;
    
    double u_bordersize;
    double u_linewidth;
    int u_gridwidth;
    double u_gridheight;
    
    long u_binBit;
    long u_logX;
    long u_logY;
    
    char *u_binAmp;
    
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
void spectro_mousemove(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers);
void spectro_mouseenter(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers);
void spectro_mouseleave(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers);
void spectro_setphase(t_spectro *x, long phase);
void spectro_realimag_perform64(t_spectro *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void spectro_dsp64(t_spectro *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
static float log_to_lin(float val, float imin, float imax, float omin, float omax);
static float lin_to_log(float val, float imin, float imax, float omin, float omax);
static float lin_to_lin(float val, float imin, float imax, float omin, float omax);
double fftAmp(int index, t_spectro *x);
double toAmp(double sample, t_spectro *x);
double todB(double sample, t_spectro *x);



static t_class *s_spectro_class;

t_symbol *ps_mousemove;
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
    class_addmethod(c, (method)spectro_mousemove,   "mousemove", A_CANT, 0); // like bidle
    class_addmethod(c, (method)spectro_dsp64,		"dsp64", A_CANT, 0);
    class_addmethod(c, (method)spectro_mouseleave,   "mouseleave",   A_CANT, 0);
    class_addmethod(c, (method)spectro_mouseenter,   "mouseenter",   A_CANT, 0);

    
    ps_mousemove = gensym("mousemove");


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
    
    CLASS_ATTR_RGBA(c, "bingrid", 0, t_spectro, u_binGrid);
    CLASS_ATTR_BASIC(c, "bingrid", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bingrid", 0, "0. 255. 0. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"bingrid",0,"rgba","Bin Grid");

	CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    /****************************************/
    
    CLASS_STICKY_ATTR(c, "category", 0, "Dimensions");
    
    CLASS_ATTR_DOUBLE(c, "bordersize", 0, t_spectro, u_bordersize);
    CLASS_ATTR_BASIC(c, "bordersize", 0);
    CLASS_ATTR_DEFAULTNAME(c, "bordersize", 0, "1.");
    CLASS_ATTR_STYLE_LABEL(c, "bordersize", 0, "double", "Border Size");
    
    CLASS_ATTR_DOUBLE(c, "linewidth", 0, t_spectro, u_linewidth);
    CLASS_ATTR_BASIC(c, "linewidth", 0);
    CLASS_ATTR_DEFAULTNAME(c, "linewidth", 0, "1.");
    CLASS_ATTR_STYLE_LABEL(c, "linewidth", 0, "double", "Line Width");
    
    CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    /****************************************/
    
    CLASS_STICKY_ATTR(c, "category", 0, "Bin-Amplitude");
    
    CLASS_ATTR_LONG(c, "binamp", 0, t_spectro, u_ampOn);
    CLASS_ATTR_BASIC(c, "binamp", 0);
    CLASS_ATTR_STYLE_LABEL(c, "binamp", 0, "onoff", "Show Bin Amplitude");
    CLASS_ATTR_SAVE(c, "binamp", 0);
    
    CLASS_ATTR_RGBA(c, "textcolor", 0, t_spectro, u_text);
    CLASS_ATTR_BASIC(c, "textcolor", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "textcolor", 0, "0. 0. 255. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"textcolor",0,"rgba","Text Color");
    
    CLASS_ATTR_RGBA(c, "linecolor", 0, t_spectro, u_verticalLine);
    CLASS_ATTR_BASIC(c, "linecolor", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "linecolor", 0, "0. 0. 255. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"linecolor",0,"rgba","Vertical Line Color");
    
    CLASS_ATTR_INT32(c, "fontsize", 0, t_spectro, u_fontSize);
    CLASS_ATTR_BASIC(c, "fontsize", 0);
    CLASS_ATTR_DEFAULTNAME(c, "fontsize", 0, "12");
    CLASS_ATTR_STYLE_LABEL(c,"fontsize",0,"int","Font Size");

    CLASS_STICKY_ATTR_CLEAR(c, "category");

    /****************************************/
    
    CLASS_STICKY_ATTR(c, "category", 0, "FFT");
    
    CLASS_ATTR_LONG(c, "bins", 0, t_spectro, u_binBit);
    CLASS_ATTR_LABEL(c, "bins", 0, "Number of Bins");
    CLASS_ATTR_DEFAULT(c, "bins", 0, "4");
    CLASS_ATTR_ENUMINDEX5(c, "bins", 0, "128", "256", "512", "1024", "2048");
    CLASS_ATTR_BASIC(c, "bins", 0);
    
    CLASS_ATTR_LONG(c, "interpX", 0, t_spectro, u_logX);
    CLASS_ATTR_LABEL(c, "interpX", 0, "Frequency Interpolation");
    CLASS_ATTR_DEFAULT(c, "interpX", 0, "0");
    CLASS_ATTR_ENUMINDEX2(c, "interpX", 0, "Linear", "Logarithmic");
    CLASS_ATTR_BASIC(c, "interpX", 0);
    
    CLASS_ATTR_LONG(c, "interpY", 0, t_spectro, u_logY);
    CLASS_ATTR_LABEL(c, "interpY", 0, "Amplitude Interpolation");
    CLASS_ATTR_DEFAULT(c, "interpY", 0, "0");
    CLASS_ATTR_ENUMINDEX2(c, "interpY", 0, "Linear", "Logarithmic");
    CLASS_ATTR_BASIC(c, "interpY", 0);
    
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
    
    //*g: background
	t_jgraphics *g = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
    //*h: samples
    t_jgraphics *h = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
    //*t freq text
    t_jgraphics *t = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
    //*l cursor line
    t_jgraphics *l = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
    //*b bin grid
    t_jgraphics *b = (t_jgraphics *) patcherview_get_jgraphics(patcherview);
    
    
    
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
    
	//paint samples
    jgraphics_set_source_jrgba(h, &x->u_samples);
    jgraphics_set_line_width(h, x->u_linewidth);
    
    
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
        
        x->u_ampOn = 0;
        
    } else {

        //number of bins below nyquist
        int numberOfBins = (int) (x->f_fftsize / 2);
        double sample[numberOfBins];
        float range;
        
        x->u_ampOn = 1;
        
        //draw spectroscope
        for(int i = 0; i < numberOfBins; i++){
            
            //scale samples to current gridwidth
            if(x->u_logX){
                range = ((float) i / numberOfBins) * (sys_getsr() / 2);
                range = log_to_lin(range, 0, (sys_getsr() / 2), 0, x->u_gridwidth);
            } else {
                range = ((float) i / numberOfBins) * x->u_gridwidth;
            }
            
            sample[i] = fftAmp(i, x);
            
            //draws amplitude range
            if(x->u_logY){
                sample[i] = todB(sample[i], x);
                //scale to gridheight
                if(sample[i] <= -70){
                    sample[i] = 0;
                } else {
                    sample[i] = lin_to_lin(sample[i], -70, 6, 0, x->u_gridheight);
                }
                
            } else {
                sample[i] = toAmp(sample[i], x) * x->u_gridheight;
            }
            
            if(i == 0){
                jgraphics_move_to(h, x->u_bordersize / 2, x->u_gridheight - (x->u_bordersize / 2) - sample[i]);
            } else {
                jgraphics_line_to(h, (float) range + (x->u_bordersize / 2), x->u_gridheight - sample[i]);
            }
        }
    }
    
    //paint samples
    jgraphics_set_line_width(h, x->u_linewidth);
    jgraphics_stroke(h);
    
    //paint bin amp text
    if(x->u_ampOn && x->u_mouseover == true){
        
        //free char array
        free(x->u_binAmp);
        
        t_bin bins;
        
        float range = (x->u_mouseX / x->u_gridwidth);
        
        //x->f_fftsize / 2 because of nyquist mirroring
        float fftResolution = (sys_getsr() / x->f_fftsize);
        
        if(x->u_logX){
            range = lin_to_log(range, 0, 1.0, 0, 1.0);
        }
        
        bins.number = range * (x->f_fftsize / 2);
        bins.freq = fftResolution * bins.number;
        bins.amp = fftAmp(bins.number, x);
        
        if(x->u_logY){
            if(bins.amp == 0){
                asprintf(&x->u_binAmp, "%iHz: -infdB", bins.freq);
            } else {
                bins.amp = todB(bins.amp, x);
                asprintf(&x->u_binAmp, "%iHz: %fdB", bins.freq, bins.amp);
            }
        } else {
            bins.amp = toAmp(bins.amp, x);
            asprintf(&x->u_binAmp, "%iHz: %f", bins.freq, bins.amp);
        }
        
        //paint text
        t_jtextlayout	*mytxt;
        t_jfont	*myfont;
        
        object_attr_getjrgba((t_object *)x, _sym_textcolor, &x->u_text);
        
        mytxt = jtextlayout_create();
        myfont = jfont_create(systemfontname(), JGRAPHICS_FONT_SLANT_NORMAL, JGRAPHICS_FONT_WEIGHT_NORMAL, 11.0);
        jfont_set_font_size(myfont, x->u_fontSize);
        jtextlayout_set(mytxt,
                        x->u_binAmp,
                        myfont,
                        -4, -(x->u_gridheight / 2) + 12,
                        rect.width,
                        rect.height,
                        JGRAPHICS_TEXT_JUSTIFICATION_RIGHT,
                        JGRAPHICS_TEXTLAYOUT_NOWRAP);
        
        jtextlayout_settextcolor(mytxt, &x->u_text);
        jtextlayout_draw(mytxt, t);
        
        jtextlayout_destroy(mytxt);
        jfont_destroy(myfont);
        
        //paint line
        jgraphics_set_source_jrgba(l, &x->u_verticalLine);
        jgraphics_set_line_width(l, 1);
        
        jgraphics_move_to(l, x->u_mouseX, x->u_bordersize / 2);
        jgraphics_line_to(l, x->u_mouseX, x->u_gridheight - (x->u_bordersize / 2));
        
        jgraphics_stroke(t);
        jgraphics_stroke(l);
        
    }

        //paint bin lines
        //jgraphics_set_source_jrgba(b, &x->u_binGrid);
        //jgraphics_set_line_width(b, 1);
        
        //jgraphics_move_to(b, /*POSITION*/, x->u_bordersize / 2);
        //jgraphics_line_to(b, /*POSITION*/, x->u_gridheight - (x->u_bordersize / 2));
            
        //jgraphics_stroke(b);

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
    
    long fftsize = FFT_MAX_BINS;
    long interval = fftsize;
    long phase = 0;
    long inverse = (s == ps_ifft);
    
    if (!(d = object_dictionaryarg(argc,argv)))
        return NULL;
    
    x = (t_spectro *)object_alloc(s_spectro_class);
    
    x->f_fftsize = fftsize; // FFT size
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

double fftAmp(int index, t_spectro *x)
{
    return (sqrt((x->f_imagout[index] * x->f_imagout[index]) + (x->f_realout[index] * x->f_realout[index])));
}

double toAmp(double sample, t_spectro *x)
{
    float mult = 1.0 / (x->f_fftsize/PI);
    if(sample > 0){
        sample *= mult;
    }
    return sample;
}

double todB(double sample, t_spectro *x)
{
    float mult = 1.0 / (x->f_fftsize/PI);
    if(sample > 0){
        sample *= mult;
        sample = (20 * log10(sample));
    }
    return sample;
}

void spectro_mouseenter(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->u_mouseover = true;
    jbox_redraw((t_jbox *)x);
}

void spectro_mouseleave(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->u_mouseover = false;
    jbox_redraw((t_jbox *)x);
}

void spectro_mousemove(t_spectro *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->u_mouseX = pt.x;
}

void spectro_setphase(t_spectro *x, long phase)
{
    x->f_phase = phase;
}

static float log_to_lin(float val, float imin, float imax, float omin, float omax)
{
    float logterm = ((val - imin) / (imax - imin)) + 0.0001f;
    if (logterm <= 0) {
        return omin;
    }
    return ((((log(logterm) / 9.21034f + 1.f) / 1.000011f) * (omax - omin)) + omin);
    
}

static float lin_to_log(float val, float imin, float imax, float omin, float omax)
{
    if ((val - imin) == 0 || (imax - imin) == 0) {
        return omin;
    }
    return (((exp(((val - imin) / (imax - imin) * 1.000011 - 1.) * 9.21034) - 0.0001) * (omax - omin)) + omin);
}


static float lin_to_lin(float val, float imin, float imax, float omin, float omax)
{
    if (val == 0 || imax == 0) {
        return 0;
    }
    return ((((val - imin) / (imax - imin)) * (omax - omin) + omin));
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


//FFT code borrowed from the fft~ object
void spectro_realimag_perform64(t_spectro *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) // based on fft_complex_perform()
{
    //check if bin argument changed
    
    //update bin attribute
    x->f_fftsize = FFT_MIN_BINS << x->u_binBit;
    
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

