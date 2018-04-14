/*
	oscope~ -- new and improved oscilloscope
 
    Part of a new group of smart audio visualization objects for Max which includes an oscilloscope,
    a spectroscope, a 3D spectroscope, smart meters, and maybe more
*/

#include "ext.h"
#include "ext_obex.h"
#include "jpatcher_api.h"
#include "jgraphics.h"
#include "z_dsp.h"  
#include <string.h>


typedef struct _oscope
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
    
    double *u_buffer;
    int u_bufferSize;
    int k_dsp;
    int k_paint;
    
    long attrLatching;
    
} t_oscope;


void *oscope_new(t_symbol *s, long argc, t_atom *argv);
void oscope_bang(t_oscope *x);
void oscope_free(t_oscope *x);
void oscope_assist(t_oscope *x, void *b, long m, long a, char *s);
void oscope_paint(t_oscope *x, t_object *patcherview);
void oscope_getdrawparams(t_oscope *x, t_object *patcherview, t_jboxdrawparams *params);
void oscope_perform64(t_oscope *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void oscope_dsp64(t_oscope *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);


static t_class *s_oscope_class;

void ext_main(void *r)
{
	t_class *c;

	c = class_new("oscope~", (method)oscope_new, (method)oscope_free, sizeof(t_oscope), 0L, A_GIMME, 0);

	c->c_flags |= CLASS_FLAG_NEWDICTIONARY;
	jbox_initclass(c, JBOX_FIXWIDTH | JBOX_COLOR);
    class_dspinitjbox(c);

	class_addmethod(c, (method)oscope_paint,		"paint",	A_CANT, 0);
	class_addmethod(c, (method)oscope_assist,		"assist",	A_CANT, 0);
    class_addmethod(c, (method)oscope_dsp64,		"dsp64", A_CANT, 0);
    class_addmethod(c, (method)oscope_bang,		    "bang", A_CANT, 0);

    
    
	// attributes
	CLASS_STICKY_ATTR(c, "category", 0, "Color");
    
	CLASS_ATTR_RGBA(c, "bgcolor", 0, t_oscope, u_background);
    CLASS_ATTR_BASIC(c, "bgcolor", 0);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bgcolor", 0, "1. 1. 1. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bgcolor",0,"rgba","Background Color");
    
    CLASS_ATTR_RGBA(c, "samplecolor", 0, t_oscope, u_samples);
    CLASS_ATTR_BASIC(c, "samplecolor", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "samplecolor", 0, "0. 0. 255. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"samplecolor",0,"rgba","Sample Color");

	CLASS_ATTR_RGBA(c, "bordercolor", 0, t_oscope, u_outline);
    CLASS_ATTR_BASIC(c, "bordercolor", 0);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bordercolor", 0, "0. 0. 0. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bordercolor",0,"rgba","Border Color");

	CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    CLASS_STICKY_ATTR(c, "category", 0, "Dimensions");
    
    CLASS_ATTR_DOUBLE(c, "bordersize", 0, t_oscope, u_bordersize);
    CLASS_ATTR_BASIC(c, "bordersize", 0);
    CLASS_ATTR_DEFAULTNAME(c, "bordersize", 0, "1.");
    CLASS_ATTR_STYLE_LABEL(c, "bordersize", 0, "double", "Border Size");
    
    CLASS_ATTR_DOUBLE(c, "linewidth", 0, t_oscope, u_linewidth);
    CLASS_ATTR_BASIC(c, "linewidth", 0);
    CLASS_ATTR_DEFAULTNAME(c, "linewidth", 0, "3.");
    CLASS_ATTR_STYLE_LABEL(c, "linewidth", 0, "double", "Line Width");
    
    CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    
    CLASS_STICKY_ATTR(c, "category", 0, "Waveform");
    
    CLASS_ATTR_LONG(c, "latching", 0, t_oscope, attrLatching);
    CLASS_ATTR_BASIC(c, "latching", 0);
    CLASS_ATTR_STYLE_LABEL(c, "latching", 0, "onoff", "Zero-Cross Latching");
    CLASS_ATTR_SAVE(c, "latching", 0);
    
    CLASS_STICKY_ATTR_CLEAR(c, "category");


	CLASS_ATTR_DEFAULT(c,"patching_rect",0, "0. 0. 400. 250.");

	class_register(CLASS_BOX, c);
	s_oscope_class = c;
}

void oscope_assist(t_oscope *x, void *b, long m, long a, char *s)
{
	if (m == 1)		//inlet
		sprintf(s, "(signal) Audio Input");
}

void oscope_paint(t_oscope *x, t_object *patcherview)
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
    
    //draw initial waveform (x = 0)
    if(sys_getdspstate() == 0){
        x->u_gridwidth = rect.width - x->u_bordersize;
        x->u_gridheight = rect.height / 2;
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
        
        x->u_gridwidth = rect.width - x->u_bordersize;
        x->u_gridheight = rect.height / 2;
        x->k_paint =((x->k_dsp - x->u_gridwidth + x->u_bufferSize) % x->u_bufferSize);
        
        //zero-cross latching
        if(x->attrLatching){
            bool isCrossing = 1;
                while(!(x->u_buffer[x->k_paint] < 0 && x->u_buffer[(x->k_paint + 1) % x->u_bufferSize] > 0)){
                    if(x->k_paint == 0){
                        x->k_paint = x->u_bufferSize;
                    }
                    x->k_paint--;
                    if(x->k_paint == x->k_dsp){
                        isCrossing = 0;
                        break;
                    }
                }
            
            double w2, w1;
            
            if(isCrossing == 0){
                w2 = 0;
                w1 = 0;
            } else {
                w2 = fabs(x->u_buffer[x->k_paint]) / fabs(x->u_buffer[x->k_paint] - x->u_buffer[(x->k_paint + 1) % x->u_bufferSize]);
                w1 = 1 - w2;
            }
        
            //draw waveform
            for(int i = 0; i < x->u_gridwidth; i++){
                int index = (x->k_paint + i) % x->u_bufferSize;
                
                double y1 = x->u_buffer[index];
                double y2 = x->u_buffer[(index + 1) % x->u_bufferSize];
                
                double sample = (y1 * w1) + (y2 * w2);
                
                //scale samples to current gridwidth
                sample = ((-sample + 1) / 2) * (rect.height - (x->u_bordersize * 2)) + x->u_bordersize;
                
                if(i == 0){
                    jgraphics_move_to(h, x->u_bordersize / 2, sample);
                } else {
                    jgraphics_line_to(h, (double) i + (x->u_bordersize / 2), sample);
                }
            }


        } else {
            
            //draw waveform
            for(int i = 0; i < x->u_gridwidth; i++){
                int index = (x->k_paint + i) % x->u_bufferSize;

                double sample = x->u_buffer[index];
                
                //scale samples to current gridwidth
                sample = ((-sample + 1) / 2) * (rect.height - (x->u_bordersize * 2)) + x->u_bordersize;
                
                if(i == 0){
                    jgraphics_move_to(h, x->u_bordersize / 2, sample);
                } else {
                    jgraphics_line_to(h, (double) i + (x->u_bordersize / 2), sample);
                }
            }

        }
    
    }
    
    jgraphics_set_line_width(h, x->u_linewidth);
    jgraphics_stroke(h);
}


void oscope_getdrawparams(t_oscope *x, t_object *patcherview, t_jboxdrawparams *params)
{
	params->d_bordercolor.alpha = 0;
	params->d_boxfillcolor.alpha = 0;
}

void oscope_free(t_oscope *x)
{
    dsp_freejbox((t_pxjbox *)x);
	jbox_free((t_jbox *)x);
}

void oscope_bang(t_oscope *x){
    jbox_redraw((t_jbox *)x);
}

void *oscope_new(t_symbol *s, long argc, t_atom *argv)
{
	t_oscope *x = NULL;
	t_dictionary *d = NULL;
	long boxflags;

	if (!(d = object_dictionaryarg(argc,argv)))
		return NULL;

	x = (t_oscope *)object_alloc(s_oscope_class);
	boxflags = 0
			   | JBOX_DRAWFIRSTIN
			   | JBOX_NODRAWBOX
			   | JBOX_DRAWINLAST
			   | JBOX_TRANSPARENT
               | JBOX_GROWBOTH
			   | JBOX_DRAWBACKGROUND
			   ;

    x->u_bufferSize = 10000;
    x->k_dsp = 0;
    x->k_paint = 0;
    x->u_buffer = malloc(sizeof(double)*x->u_bufferSize);
    
	jbox_new((t_jbox *)x, boxflags, argc, argv);
    x->u_box.z_box.b_firstin = (void *)x;
    dsp_setupjbox((t_pxjbox *)x,1);
    outlet_new((t_object *)x, "signal");
	attr_dictionary_process(x,d);
	jbox_ready((t_jbox *)x);
	return x;
}

void oscope_perform64(t_oscope *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    
    double    *in = ins[0];     // first inlet
    double    *out = outs[0];   // first outlet
    double     n = sampleframes; // vector size
    
    while (n--) {               // perform calculation on all samples
        if(x->k_dsp == x->u_bufferSize){
            x->u_buffer[x->k_dsp] = *in++;
            *out++ = x->u_buffer[x->k_dsp];
            x->k_dsp = 0;
        } else {
            x->u_buffer[x->k_dsp] = *in++;
            *out++ = x->u_buffer[x->k_dsp];
            (x->k_dsp)++;
        }
    }

        jbox_redraw((t_jbox *)x);

}


void oscope_dsp64(t_oscope *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    post("my sample rate is: %f", samplerate);
    
    object_method(dsp64, gensym("dsp_add64"), x, oscope_perform64, 0, NULL);

}

