/*
	oscope~ -- new and improved oscilloscope
 
 TODO:
 
 -Fix first sample's slope
 -Adjust bordersize
 -Get exact number of pixels in monitor for buffer size
 -Fixing latency of waveform (periodicity analysis, sampling at constant rate)
 -Clean up code, finalize attributes, maybe give it a grid
 -Start spectroscope
 

*/

#include "ext.h"
#include "ext_obex.h"						// required for new style Max object
#include "jpatcher_api.h"
#include "jgraphics.h"
#include "z_dsp.h"  

typedef struct _oscope
{
	t_pxjbox u_box;						// header for UI objects
	void *u_out;						// outlet pointer
	t_jrgba u_outline;
	t_jrgba u_background;
    t_jrgba u_samples;
    
    double u_boardersize;
    double u_linewidth;
    int u_gridwidth;
    double u_gridheight;
    
    double *u_buffer;
    int u_bufferSize;
    int k_dsp;
    int k_paint;
    
    bool init;
} t_oscope;


void *oscope_new(t_symbol *s, long argc, t_atom *argv);
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

    
    
	// attributes
	CLASS_STICKY_ATTR(c, "category", 0, "Color");
    
	CLASS_ATTR_RGBA(c, "bgcolor", 0, t_oscope, u_background);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bgcolor", 0, "1. 1. 1. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bgcolor",0,"rgba","Background Color");
    
    //TODO: Sample color not showing up?
    CLASS_ATTR_RGBA(c, "samplecolor", 0, t_oscope, u_samples);
    CLASS_ATTR_BASIC(c, "samplecolor", 0);
    CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "samplecolor", 0, "0. 255. 0. 1.");
    CLASS_ATTR_STYLE_LABEL(c,"samplecolor",0,"rgba","Sample Color");

	CLASS_ATTR_RGBA(c, "bordercolor", 0, t_oscope, u_outline);
	CLASS_ATTR_DEFAULTNAME_SAVE_PAINT(c, "bordercolor", 0, "0. 0. 0. 1.");
	CLASS_ATTR_STYLE_LABEL(c,"bordercolor",0,"rgba","Border Color");

	CLASS_STICKY_ATTR_CLEAR(c, "category");
    
    CLASS_STICKY_ATTR(c, "category", 0, "Dimensions");
    
    CLASS_ATTR_DOUBLE(c, "boardersize", 0, t_oscope, u_boardersize);
    CLASS_ATTR_DEFAULTNAME(c, "boardersize", 0, "1.");
    CLASS_ATTR_STYLE_LABEL(c, "boardersize", 0, "double", "Boarder Size");
    
    CLASS_ATTR_DOUBLE(c, "linewidth", 0, t_oscope, u_linewidth);
    CLASS_ATTR_DEFAULTNAME(c, "linewidth", 0, "3.");
    CLASS_ATTR_STYLE_LABEL(c, "linewidth", 0, "double", "Line Width");
    
    CLASS_STICKY_ATTR_CLEAR(c, "category");

	CLASS_ATTR_DEFAULT(c,"patching_rect",0, "0. 0. 20. 20.");

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
	jbox_get_rect_for_view((t_object *)x, patcherview, &rect);


	// paint border
	jgraphics_set_source_jrgba(g, &x->u_outline);
	jgraphics_set_line_width(g, x->u_boardersize);
	jgraphics_rectangle(g, 0. + x->u_boardersize / 2,
                        0. + x->u_boardersize / 2,
                        rect.width - x->u_boardersize,
                        rect.height - x->u_boardersize);
    
	jgraphics_stroke(g);
    
    //paint background
    jgraphics_set_source_jrgba(g, &x->u_background);
    jgraphics_rectangle(g, 0. + x->u_boardersize / 2,
                        0. + x->u_boardersize / 2,
                        rect.width - x->u_boardersize,
                        rect.height - x->u_boardersize);
    
    jgraphics_fill(g);
    
    
	// paint samples
    jgraphics_set_source_jrgba(g, &x->u_samples);
    
    //draw initial waveform (x = 0)
    if(x->init == 0){
        x->init = 1;
        x->u_gridwidth = rect.width - x->u_boardersize;
        x->u_gridheight = rect.height / 2;
        x->u_buffer = malloc(sizeof(double)*x->u_bufferSize);
        for(int i = 0; i < x->u_bufferSize; i++){
            jgraphics_rectangle(g,
                                (double) i / 10 + (x->u_boardersize / 2),
                                x->u_gridheight,
                                0.1,
                                x->u_linewidth);
        }
    } else {
        //TODO: Fix first sample's slope
        if(x->u_gridwidth != rect.width - x->u_boardersize){
            x->u_gridwidth = rect.width - x->u_boardersize;
        }
        if(x->u_gridheight != rect.height){
            x->u_gridheight = rect.height / 2;
        }
        //jgraphics_move_to(g, x->u_boardersize / 2, x->u_gridheight);
        x->k_paint = x->k_dsp;
        x->u_buffer[x->k_paint] = ((x->u_buffer[x->k_paint] + 1) / 2) * rect.height;
        
        while(!(x->u_buffer[x->k_paint] < 0 && x->u_buffer[x->k_paint + 1] > 0)){
            if(x->k_paint == 0){
                x->k_paint = x->u_bufferSize - 1;
            }
            x->k_paint--;
        }
        
        //draw waveform
        for(int i = 0; i < x->u_gridwidth; i++){
            if(x->k_paint == 0){
                x->k_paint = x->u_bufferSize - 1;
            }
            x->u_buffer[x->k_paint] = ((x->u_buffer[x->k_paint] + 1) / 2) * rect.height;
            if(i == 0){
                jgraphics_move_to(g, x->u_boardersize / 2, x->u_buffer[x->k_paint]);
                x->k_paint--;
            } else {
                jgraphics_line_to(g, (double) i + (x->u_boardersize / 2), x->u_buffer[x->k_paint]);
                x->k_paint--;
            }
        }
            
    }
    
    jgraphics_stroke(g);
    
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
			   //		| JBOX_NOGROW
			   //       | JBOX_GROWY
               | JBOX_GROWBOTH
			   //		| JBOX_HILITE
			   //		| JBOX_BACKGROUND
			   | JBOX_DRAWBACKGROUND
			   //		| JBOX_NOFLOATINSPECTOR
			   //		| JBOX_TEXTFIELD
			   //		| JBOX_MOUSEDRAGDELTA
			   //		| JBOX_TEXTFIELD
			   ;

    //TODO: exact buffer size of # of pixels in monitor?
    x->u_bufferSize = 10000;
    x->k_dsp = 0;
    x->k_paint = 0;
    
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
    
    if(x->init ==0){
        x->init = 1;
    }
    
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
    
    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling function
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method
    
    object_method(dsp64, gensym("dsp_add64"), x, oscope_perform64, 0, NULL);

}

