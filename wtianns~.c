//
//  wtianns~.c
//  wtianns~
//  wavetable-inspired artificial neural network synthesis
//  this version implements miller puckette's phase-bashed packet synthesis.
//  
//  Created by Nicholas Solem on 9/1/21.
//  Copyright © 2021 Nicholas Solem. All rights reserved.
//
#include "m_pd.h"
#define WINDOWING 0
#define _180ON 1
#define OUTPUT_SIZE 40


//#define K2C
#ifndef K2C
#include "fann.h"
#else
#include "k2c_tensor_include.h"
#define M_PI 3.141592653589793 // delete once i link <model>.o (?) because <model>.c includes <math.h>
#endif

#include <fftw3.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

//#define WAVELENGTH 1024

t_class *wtianns_tilde_class;

typedef struct wtianns_tilde
{
    t_object x_obj;     /* obligatory header */
    t_float intern_freq;
    
#ifndef K2C
    struct fann *ann;
    fann_type *calc_out;
#else
    float *calc_out;
#endif
    
    int winSize;
    double phase;
    float fs, fs_inv; // inverse of sample rate
    float T, S;     // params: T controls windowing function width (corresponding to bandwidth); S controls phasor multiplier (corresponding with hard sync)
    
    // frequency-domain wavetable, DC to Nyquist only, phase-bashed
    //double *freqtable;//
    //double *wavetable;// [WAVELENGTH];
    
    //  a temporary bool that set on every 'knob' change & resets at end of block
    int tabchanged;
    // a bool that flips at every 'knob' change
    int whichTab;
    
    float *freqtable, *scratchTable;
    float *bin_mins, *bin_ranges; // denormalization
    // alternates writing to A and B so we can xfade
    float *wavetable_A;
    float *wavetable_B;
    float *winFunc;
    fftwf_plan thePlan;
    
    t_float testthing;
    
    t_outlet *main_out;
    t_outlet *out2, *out3, *out4, *out5;
    
} t_wtianns_tilde;
//NEED TO APPLY TWICE TO END UP WITH 4W WORTH
// IDCT GIVES HALF-WAVE; PHASEBASH METHOD REQUIRES 2 FULL PERIODS PRIOR TO WINDOWING
// if freq of original if 1/W, then Hann window will give freqs of 1/2W and 3/2W (aliased frequencies)
// those particular freqs will get cancelled out when you overlap-add them
// whereas before there was no cancellation when using single period
// could have instead insured that we are only writing every other bin
// by default, f0 should occupy bin 1, then user can change to have subharmonics if they choose
void xt_halfspec2fullspec(int winSize, float *halfspec, float *fullspec) {
    int _2W = winSize * 2;
    for (int n = 0; n < _2W; ++n) {
        int inIdx;
        if (n < winSize) // [0..2047]
            inIdx = n;
        else            // [2048..4095]
            inIdx = _2W - n - 1;
        fullspec[n] = halfspec[inIdx];
    }
}
void squarespec(int winSize, float *inSpec, float *outSpec) {
    for (int n=0; n < winSize; n++){
        outSpec[n] = inSpec[n] * inSpec[n];
    }
}
float maxVal(float *arr, float *weightTab, int N)
{
    float max = 0.f;
    float tmp;
    // STARTS INDEXING AT 2ND VAL, ONLY NEEDS TO GO TO HALFPOINT OF WAVEFORM
    int n,m;
    for (n = 1; n < N/2; n++)
    {
        m = (n + N/2) % N;
        tmp = arr[n] * weightTab[n] + arr[m] * weightTab[m];
        tmp = tmp > 0.f ? tmp : -tmp;
        if (tmp > max)
            max = tmp;
    }
    return max;
}
float cubicInterp(float *wavetable, float fIdx, int winSize)
{// adapted from https://www.musicdsp.org/en/latest/Other/49-cubic-interpollation.html?highlight=cubic
    int iIdx = (int)fIdx;// assuming we never get a negative fIdx; then it would round up
    float frac = fIdx - (float)iIdx;
    
    float ym1 = wavetable[(iIdx - 1) & (winSize - 1)];    // could &= the mask here
    float y0  = wavetable[(iIdx + 0) & (winSize - 1)];
    float y1  = wavetable[(iIdx + 1) & (winSize - 1)];    // and here...
    float y2  = wavetable[(iIdx + 2) & (winSize - 1)];    // and here...
    float a = (3.f * (y0-y1) - ym1 + y2) * 0.5f;
    float b = 2.f*y1 + ym1 - (5.f*y0 + y2) * 0.5f;
    float c = (y1 - ym1) * 0.5f;
    float y = (((a * frac) + b) * frac + c) * frac + y0;
    return y;
}
void wtianns_tilde_calc_net(t_wtianns_tilde *x, t_symbol *s, int argc, t_atom *argv) {
    if (x->ann == NULL){
        error("ann not initialized");
        return; }
    unsigned int inputSize;
#ifndef K2C
    inputSize = x->ann->num_input;
#else
#endif
    if(argc < inputSize) {
        error("wtianns~: too few input values!!");
        return; }
    float knobs[inputSize];
    
    int i;
    float tmp;
    for (i = 0; i < inputSize; i++) {
        tmp = atom_getfloat(argv++);
        knobs[i] = tmp;
    }
    //int zeroOutStartBin = (int)((1.f / (knobs[0] * 2.f)) * x->winSize);
#ifndef K2C
    x->calc_out = fann_run(x->ann, (fann_type *)knobs);
#else
    x->calc_out = 
#endif
    //xt_halfspec2fullspec(x->winSize, x->calc_out, x->freqtable);

    // GOOD TRICK TO TRY BASED ON TRIANGLE WAVE:
    // INVERT EVERY OTHER HARMONIC SO SUMMATION AT 0 PHASE IS MINIMIZED
    // could also try rotating phase as [0, pi/2, pi, 3/2pi] if using DFT
    // try inverting every 4th harmonic to get rid of peak
    for (i=0; i < OUTPUT_SIZE; i++){
        x->freqtable[i] = (x->calc_out[i] * x->bin_ranges[i]) + x->bin_mins[i];//0.f;
    }
//    for (i=zeroOutStartBin; i < x->winSize; i++)
//        x->freqtable[i] = 0.f;
    squarespec(x->winSize, x->freqtable, x->freqtable);
    //NEED TO DO THIS TWICE
    squarespec(x->winSize, x->freqtable, x->freqtable);
    x->freqtable[0] = 0.f; // block DC
    //x->freqtable[1] = 1.f;
//    x->freqtable[1] = x->calc_out[1];// 1.f;
//    x->freqtable[2] = x->calc_out[2];// 1.f;
//    x->freqtable[3] = x->calc_out[3];// 1.f;
    //x->freqtable[3] = 1.f;
    x->tabchanged = 1;

    float maxv;
    fftwf_execute(x->thePlan);// fft from freqtable to wavetable
    float *targetTable;
    if (x->whichTab == 1) {// time 0, we want to WRITE to B while it's been READING from A
        targetTable = x->wavetable_B;
    }
    else {
        targetTable = x->wavetable_A;
    }

    xt_halfspec2fullspec(x->winSize, x->scratchTable, targetTable);     // reflect half-cycle to single cycle
    xt_halfspec2fullspec(x->winSize * 2, targetTable, targetTable);     // reflect single cycle to double cycle

    maxv = maxVal(targetTable, x->winFunc, x->winSize * 4);       // now enough info to find max (windowing function peaks )
    x->testthing = maxv;

    for (i = 0; i < x->winSize * 4; i++){
        targetTable[i] = targetTable[i] / maxv;
    }
    
    // alternate which to read from for _perform method
    x->whichTab += 1;
    x->whichTab %= 2;
}

static t_int *wtianns_tilde_perform(t_int *w)
{
    t_wtianns_tilde *x =  (t_wtianns_tilde *)(w[1]);
    t_float *in =       (t_float *)(w[2]); // NEED TO INCREMENT POINTER EACH SAMPLE
    t_float *out =      (t_float *)(w[3]);
    t_float *out2 =      (t_float *)(w[4]);
    t_float *out3 =      (t_float *)(w[5]);
    t_float *out4 =      (t_float *)(w[6]);
    t_float *out5 =      (t_float *)(w[7]);
    int bufferlength =  (int)(w[8]);
    
    // since the newly calculated table lives as a spectrum in x->freqtable,
    // we only need to do:
    //      (at beginning of block):
    //          -check if new IFFT is needed (did 'knobs' change?) thereby set '_tabchanged' condition
    //          -cue next table by writing IFFT to it and set reference targetTable
    //      (over course of block):
    //          -if '_tabchanged', crossfade over block to target table
    //          -else continue reading from current table
// BEGINNING OF BLOCK
    
// LOOK INTO HERMITE FILTER FOR OVERSAMPLING
// hermite would interpolate on top of the oversampling
// msp uses lagrange interpolation instead of hermite, which is lower distortion under certain circumstances
// hermite has better bass frequency response than others
// lagrange may have other unpleasant artifacts, tom says
    if (x->calc_out != NULL)
    {
        // oldest, newest written tables
        float *wtLast = NULL, *wtNext = NULL;
        //
        float y0, y180, y0cubic, y180cubic;
        float win0, win180;
        float wIdx_f;

        int wIdx0, wIdx180;

        if (x->whichTab == 1){
            wtLast = x->wavetable_B;
            wtNext = x->wavetable_A;
        }
        else {
            wtLast = x->wavetable_A;
            wtNext = x->wavetable_B;
        }
        /* if (x->wtLast == NULL) x->wtLast = x->wavetable_A;
        if (x->wtNext == NULL) x->wtNext = x->wavetable_B;
        */
        int _4W = 4 * x->winSize;
        for (int buffIdx = 0; buffIdx < bufferlength; buffIdx++)
        {
            float freq = *(in++) * 0.5f;
            freq = freq < x->fs * 0.5f ? freq : x->fs * 0.5f; // highest f0 = Nyquist
            
            x->phase += freq * x->fs_inv;                // 'raw' phasor 0-1
            wIdx_f = x->phase * (float)_4W;// (x->winSize);//
            wIdx_f *= x->S * 2.f;
            
            float wIdx_f180 = wIdx_f + (float)(_4W * 0.5);//(x->winSize * 0.5);
            if (wIdx_f180 > (float)(_4W - 1))//(x->winSize - 1))
                wIdx_f180 -= (float)(_4W);//(x->winSize);

            win0 = cubicInterp(x->winFunc, wIdx_f,_4W);
            win180 = cubicInterp(x->winFunc, wIdx_f180,_4W);

            wIdx0 = (int)wIdx_f;    // floor (int) part
            wIdx0 &= _4W - 1;// x->winSize - 1; // newly added, no need for if()

            wIdx180 = wIdx0 + (int)(_4W * 0.5);//(x->winSize * 0.5);
            //if (wIdx180 >= x->winSize) wIdx180 &= x->winSize - 1;
            wIdx180 &= _4W - 1;
            float xfade = 0.f;
            if (x->tabchanged)
            { // crossfade over block
                xfade = (float)(buffIdx) / (float)(bufferlength);
                y0 =    (1.f-xfade)*wtLast[wIdx0]   + xfade*wtNext[wIdx0];
                y180 =  (1.f-xfade)*wtLast[wIdx180] + xfade*wtNext[wIdx180];
                
                y0cubic = (1.f-xfade)*cubicInterp(wtLast, wIdx_f,_4W) + (xfade)*cubicInterp(wtNext, wIdx_f,_4W);
                //if (WINDOWING) y0cubic *= win0;
                y180cubic =(1.f-xfade)*cubicInterp(wtLast, wIdx_f180,_4W) + (xfade)*cubicInterp(wtNext, wIdx_f180,_4W);
                //if (WINDOWING) y180cubic *= win180;
            }
                else // no crossfade
            {
                y0 = wtNext[wIdx0];
                y180 = wtNext[wIdx180];

                y0cubic = cubicInterp(wtNext, wIdx_f,_4W);
                //if (WINDOWING) y0cubic *= win0;
                y180cubic = cubicInterp(wtNext, wIdx_f180,_4W);
                //if (WINDOWING) y180cubic *= win180;
            }
            *out++ = (t_float)(y0cubic * win0 + y180cubic * win180);
            *out2++ = (t_float)(y0cubic);
            *out3++ = (t_float)(y180cubic);
            *out4++ = (t_float)(win0);
            *out5++ = (t_float)(x->testthing);
            while (x->phase > 1.f)
                x->phase -= 1.f;
            while (x->phase < 0.f)
                x->phase += 1.f;
        }// end of block for loop
        if (x->tabchanged == 1) x->tabchanged = 0; // MAKE SURE THAT THIS WONT INTERFERE WITH CHANGING TABLE TOWARD END OF BLOCK
    }
    else // x->calc_out (raw half spectrum) is NULL
        while (bufferlength--){
            *out++ = 0.f;
            *out2++= 1.f;
            *out3++= 1.f;
            *out4++= 1.f;
            *out5++= 1.f;}
    return (w + 9);
}

static void wtianns_tilde_dsp(t_wtianns_tilde *x, t_signal **sp) {
    x->fs = sp[0]->s_sr;
    x->fs_inv = 1.0f / sp[0]->s_sr;
    
    dsp_add(wtianns_tilde_perform,
            8,              // num items
            x,
            sp[0]->s_vec,   // input vector
            sp[1]->s_vec,   // output_vector
            sp[2]->s_vec,   // output_vector 2
            sp[3]->s_vec,   // output_vector 2
            sp[4]->s_vec,   // output_vector 2
            sp[5]->s_vec,   // output_vector 2
            sp[0]->s_n);    // blocksize
}

static void *wtianns_tilde_new(t_floatarg winSize)
{
    t_wtianns_tilde *x = (t_wtianns_tilde *)pd_new(wtianns_tilde_class);
    x->main_out = outlet_new(&x->x_obj, gensym("signal"));
    x->out2 = outlet_new(&x->x_obj, gensym("signal"));
    x->out3 = outlet_new(&x->x_obj, gensym("signal"));
    x->out4 = outlet_new(&x->x_obj, gensym("signal"));
    x->out5 = outlet_new(&x->x_obj, gensym("signal"));

    x->winSize = (int)winSize;
    //x->winSize = WAVELENGTH;
    x->whichTab = 0;
    //x->freqtable = malloc(sizeof(float *) * winSize);
    //x->wavetable = malloc(sizeof(float *) * winSize);
    float tmp;

    int _4W = x->winSize * 4;
    x->freqtable = (float *)malloc(sizeof(float) * x->winSize);
    x->bin_mins = (float *)malloc(sizeof(float) * x->winSize + 1); // a table of constant weights to denormalize the spectra
    x->bin_ranges = (float *)malloc(sizeof(float) * x->winSize + 1); // a table of constant weights to denormalize the spectra
    x->scratchTable = (float *)malloc(sizeof(float) * x->winSize); // IDCT to time domain, need to reflect waveform
    // alternates writing to A and B so we can xfade
    x->wavetable_A = (float *)malloc(sizeof(float) * _4W);
    x->wavetable_B = (float *)malloc(sizeof(float) * _4W);
    x->winFunc = (float *)malloc(sizeof(float) * _4W);
    
    x->S = 1.f;
    x->T = 1.f;
    int i;
    for (i = 0; i < _4W; i++) {
        tmp = sinf((M_PI * i) / (float)_4W);
        x->winFunc[i] = tmp * tmp;
    }
    for (i = 0; i < x->winSize+1; i++){
        x->bin_mins[i] = 0.f;
        x->bin_ranges[i] = 1.f;
    }
    
    //x->thePlan = fftwf_plan_r2r_1d(x->winSize, x->freqtable, x->scratchTable, FFTW_BACKWARD, FFTW_MEASURE);
    x->thePlan = fftwf_plan_r2r_1d(x->winSize, x->freqtable, x->scratchTable,  FFTW_REDFT01, FFTW_MEASURE);
#ifdef K2C
    deletethismodel_initialize();
#endif
    x->calc_out = NULL;
    x->ann = NULL;
    
    return (x);
}
static void wtianns_tilde_set_S(t_wtianns_tilde *x, float S){
    x->S = S;
}
static void wtianns_tilde_set_T(t_wtianns_tilde *x, float T){
    x->T = T;
}
    // fill this out and attach as destructor
    static void wtianns_tilde_delete(t_wtianns_tilde *x) {
        if (x->ann != NULL)
            fann_destroy(x->ann);
        free((void *)x->freqtable);
        free((void *)x->bin_mins);
        free((void *)x->bin_ranges);
        free((void *)x->scratchTable);
        free((void *)x->wavetable_A);
        free((void *)x->wavetable_B);
        free((void *)x->winFunc);
#ifdef K2C
        deletethismodel_terminate();
#endif
    }
#ifndef K2C
    static void wtianns_tilde_read(t_wtianns_tilde *x, t_symbol *s, int argc, t_atom *argv) {
        const char *stringy = argv->a_w.w_symbol->s_name;
        if (x->ann != NULL)
            fann_destroy(x->ann);
        
        x->ann = fann_create_from_file(stringy);
        if (x->ann != NULL) {
            const unsigned int numNeurons = fann_get_total_neurons(x->ann);
            post("num neurons: %i", numNeurons);
        }
        else
            post("read from file unsuccessful");
    }
    static void wtianns_tilde_load_bin_weights(t_wtianns_tilde *x, t_symbol *s, int argc, t_atom *argv) { // .coef file
        const char *fn = argv->a_w.w_symbol->s_name;
        FILE *bwf = fopen(fn, "r");
        
        if (bwf != NULL) {
            char str[1024];
            /* assumes no word exceeds length of 1023 */
            float line_d;
            int i = 0;
            while (fscanf(bwf, " %1023s", str) == 1) {
                line_d = atof(str);
                if (i < x->winSize+1)
                    x->bin_mins[i] = line_d;
                else
                    x->bin_ranges[i - x->winSize - 1] = line_d;
                //post("float %f", (float)(line_d));
                i++;
            }
        }
        else
            post("file NULL");
        int status = fclose(bwf);
        status = status;
    }
#endif
    
    void wtianns_tilde_write(t_wtianns_tilde *x){
        FILE *fp;
        fp = fopen("/Users/nicholassolem/development/wtianns~/wtiannsTest.txt", "w");
        if (x->calc_out != NULL)
        {
            for (int i = 0; i < x->winSize; i++) {
                //fprintf(fp, "%i, \t%f,\t%f,\t%f\n", i, x->freqtable[i], x->wavetable_A[i], x->wavetable_B[i]);
                fprintf(fp, "%f, ", x->freqtable[i]);
            }
        }
        fclose(fp);
    }
    void wtianns_tilde_setup(void)
    {
        wtianns_tilde_class = class_new(gensym("wtianns~"), (t_newmethod)wtianns_tilde_new, (t_method)wtianns_tilde_delete, sizeof(t_wtianns_tilde), 0, A_DEFFLOAT, 0);
#ifndef K2C
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_read, gensym("read"), A_GIMME, 0);
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_load_bin_weights, gensym("load_bin_weights"),A_GIMME,0);
#endif
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_set_S, gensym("S"), A_DEFFLOAT, 0);
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_set_T, gensym("T"), A_DEFFLOAT, 0);
        //class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_calc_net, gensym("calc"), 0);
        class_addlist(wtianns_tilde_class, (t_method)wtianns_tilde_calc_net);
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_write, gensym("write"), 0);
        /* this is magic to declare that the leftmost, "main" inlet
         takes signals; other signal inlets are done differently... */
        CLASS_MAINSIGNALIN(wtianns_tilde_class, t_wtianns_tilde, intern_freq);
        class_addmethod(wtianns_tilde_class, (t_method)wtianns_tilde_dsp, gensym("dsp"), 0);
        class_sethelpsymbol(wtianns_tilde_class, gensym("wtianns~"));
    }

