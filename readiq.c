/*
 * readiq.c
 *
 * C-MEX function for MATLAB to read IQ data from PX-1000
 * This is specific to certain files only. PX-1000's software 
 * is still changing
 *
 * Copyright 2012 Boon Leng Cheong. All rights reserved.
 *
 */
#include <math.h>
#include "mex.h"

#include "net.h"
#include "iqpulse.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* Math macros */
#ifndef NAN
#define NAN          (0.0/0.0)
#endif
#ifndef MAX
#define MAX(X,Y)     ((X)>(Y)?(X):(Y))
#endif
#ifndef MIN
#define MIN(X,Y)     ((X)>(Y)?(Y):(X))
#endif
#ifndef SQ
#define SQ(X)        ((X)*(X))
#endif
#ifndef ABS
#define ABS(X)       ((X)<0?-(X):(X))
#endif

/**************************************************
 *
 *  i n t 2 s t r 3
 *
 **************************************************
 *
 *  Integer to string with 3-digit grouping
 */
static void int2str3(char *linebuff, long long num) {
	int idx,jdx,kdx;
	
	sprintf(linebuff,"%lld",num);
	if (num<1000)
		return;
	else {
		/*kdx = (int)floor(log10((double)num)/3);*/
		kdx = (int)(strlen(linebuff)-1)/3;
		idx = (int)strlen(linebuff)+kdx; jdx = 1;
		linebuff[idx] = '\0';
		while (idx>0) {
			idx--;
			linebuff[idx] = linebuff[idx-kdx];
			if (jdx>3) {
				jdx = 0;
				linebuff[idx] = ',';
				kdx--;
			}
			jdx++;
		}
	}
}



/**************************************************
 *
 *  d a t e s t r
 *
 **************************************************/
static void datestr(char *timestr, const struct timeval *tnum, const int maxchar, const int style) {
	switch (style) {
		case 0:
			/* 2006/10/03 11:30:00.123456 */
			strftime(timestr,(size_t)maxchar-8,"%Y/%m/%d %H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			sprintf(timestr+strlen(timestr),".%06d",(int)(tnum->tv_usec));
			break;
		case 1:
			/* 2006/10/03 11:30:00.123 */
			strftime(timestr,(size_t)maxchar-5,"%Y/%m/%d %H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			sprintf(timestr+strlen(timestr),".%03d",(int)(tnum->tv_usec/1000));
			break;
		case 2:
			/* 11:30:00.123456 */
			strftime(timestr,(size_t)maxchar-8,"%H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			sprintf(timestr+strlen(timestr),".%06d",(int)tnum->tv_usec);
			break;
		case 3:
			/* 11:30:00.123 */
			strftime(timestr,(size_t)maxchar-5,"%H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			sprintf(timestr+strlen(timestr),".%06d",(int)tnum->tv_usec);
		case 4:
			/* 2006/10/03 11:30:00 */
			strftime(timestr,(size_t)maxchar,"%Y/%m/%d %H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			break;
		case 5:
			/* 2006/10/03 */
			strftime(timestr,(size_t)maxchar,"%Y/%m/%d",gmtime((time_t*)&(tnum->tv_sec)));
			break;
		case 6:
			/* 11:30:00 */
			strftime(timestr,(size_t)maxchar,"%H:%M:%S",gmtime((time_t*)&(tnum->tv_sec)));
			break;
		case 7:
			/* YYYYMMDD-HHMMSS */
			strftime(timestr,(size_t)maxchar,"%Y%m%d-%H%M%S",gmtime((time_t*)&(tnum->tv_sec)));
			break;
		default:
			timestr[0] = '\0';
			break;
	}
	return;
}

// Global
int
	run = 1;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	FILE                 *fid;
	char                 linebuff[1000], filename[160], msg[2048];
	unsigned char        *iqbuff;
	size_t               fsize, read_size, rx_byte;
	int                  i, j, k, l, read_count, ip;

	struct packet_header ph;

	mxArray              *tmpmx, *iq_struct, *ch1, *ch2, *elmx, *azmx;
	int                  num_pulses = 1000;
	float                *ch1_i, *ch1_q, *ch2_i, *ch2_q;

	int                  c;
	int                  pulse_size = 0;
	PXIQPulse            pulse;
	PXIQFileHeader       file_header;

#define MAX_DWELL  256

	if (nrhs < 1)
		mexErrMsgTxt("I need one input (vector or matrix)."); 
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Filename must be a string.");

	mxGetString(prhs[0], filename, 159);

#define MAX_PULSES_PER_FILE 60000

	if (nrhs > 1) {
		k = (int)mxGetScalar(prhs[1]);
		mexPrintf("Forward to pulse %d\n", k);
	} else {
		k = 0;
	}
	if (nrhs > 2) {
		l = (int)mxGetScalar(prhs[2]);
		mexPrintf("Read %d pulses\n", l);
	} else {
		l = MAX_PULSES_PER_FILE;
	}

	// IQ of one radial
	iqbuff = (unsigned char *)mxMalloc(PX_MAX_GATES*sizeof(int16_t)*2);

	if ((fid = fopen(filename, "rb")) == NULL)
		mexErrMsgTxt("I cannot open the file %s\n");
	fseek(fid, 0, SEEK_END);
	fsize = ftell(fid);
	rewind(fid);
	
	fread(&file_header, sizeof(PXIQFileHeader), 1, fid);
    printf("header.build = %d\n", file_header.build);
    printf("header.radar = %s\n", file_header.radar);
    printf("header.radar_lat_deg = %f\n", file_header.radar_lat_deg);
    printf("header.radar_lon_deg = %f\n", file_header.radar_lon_deg);
	
    // Read the very first pulse's header
    fseek(fid, sizeof(struct packet_header), SEEK_CUR);
	fread(&pulse, sizeof(PXIQPulseHeader), 1, fid);
	printf("ngate:%d  pw:%d\n", (int)pulse.ngate, (int)pulse.pw_n);

    // A fixed ngate for the rest, not a good idea
	const int ngate = pulse.ngate;

 	int packet_size = sizeof(struct packet_header) + sizeof(PXIQPulseHeader) + 2*ngate*sizeof(PXIQ);
 	num_pulses = (fsize - sizeof(PXIQFileHeader))/packet_size;
    printf("packet_size:%d\n", packet_size);
 	printf("Estimated to have %d pulse(s). (%.2f)\n", num_pulses, (float)(fsize - sizeof(PXIQFileHeader))/packet_size);

    // Reverse back to end of file header
	fseek(fid, sizeof(PXIQFileHeader), SEEK_SET);
	printf("%d vs %d\n", ftell(fid), sizeof(PXIQFileHeader));
	
	if (num_pulses > MAX_PULSES_PER_FILE) {
		mexPrintf("Exceed maximum of pulses I can handle. Truncate to %d\n", l);
	}
   	num_pulses = MIN(num_pulses, l);

//	fseek(fid, k*(sizeof(PXIQPulse) + sizeof(struct packet_header)), SEEK_SET);

	int2str3(linebuff, fsize);
	mexPrintf("Filename: %s ( %s B ) %d %d\n", filename, linebuff, packet_size, fsize-sizeof(PXIQFileHeader));
#ifdef DEBUG
	mexPrintf("sizeof(struct packet_header) = %d\n", (int)sizeof(struct packet_header));
#endif
	const char *iqfields[] = {"ch1", "ch2","radar","task","waveform","lat_deg","lon_deg","delr_m",
        "start_gate","filter_size1","filter_size2","el_deg","az_deg","el_int","az_int"};
	iq_struct = mxCreateStructMatrix(1, 1, sizeof(iqfields)/sizeof(iqfields[0]), (const char**)iqfields);
	plhs[0] = iq_struct;

    // Now, we populate the structure

    ch1 = mxCreateNumericMatrix(ngate, num_pulses, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetField(iq_struct, 0, "ch1", ch1);
    ch1_i = (float *)mxGetData(ch1);
    ch1_q = (float *)mxGetImagData(ch1);

    ch2 = mxCreateNumericMatrix(ngate, num_pulses, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetField(iq_struct, 0, "ch2", ch2);
    ch2_i = (float *)mxGetData(ch2);
    ch2_q = (float *)mxGetImagData(ch2);
					
	azmx = mxCreateNumericMatrix(1, num_pulses, mxSINGLE_CLASS, mxREAL);
	float *az = (float *)mxGetData(azmx);
	
	elmx = mxCreateNumericMatrix(1, num_pulses, mxSINGLE_CLASS, mxREAL);
	float *el = (float *)mxGetData(elmx);
    
	mxSetField(iq_struct, 0, "radar", mxCreateString(file_header.radar));
    mxSetField(iq_struct, 0, "task", mxCreateString(file_header.task));
    mxSetField(iq_struct, 0, "waveform", mxCreateString(file_header.waveform));
	mxSetField(iq_struct, 0, "lat_deg", mxCreateDoubleScalar(file_header.radar_lat_deg));
	mxSetField(iq_struct, 0, "lon_deg", mxCreateDoubleScalar(file_header.radar_lon_deg));
    mxSetField(iq_struct, 0, "delr_m", mxCreateDoubleScalar(30.0));
    mxSetField(iq_struct, 0, "start_gate", mxCreateDoubleScalar((double)file_header.start_gate+1));
    mxSetField(iq_struct, 0, "filter_size1", mxCreateDoubleScalar((double)file_header.filter_size1));
    mxSetField(iq_struct, 0, "filter_size2", mxCreateDoubleScalar((double)file_header.filter_size2));
    mxSetField(iq_struct, 0, "el_deg", elmx);
    mxSetField(iq_struct, 0, "az_deg", azmx);
    
    tmpmx = mxCreateNumericMatrix(1, num_pulses, mxINT16_CLASS, mxREAL);
    mxSetField(iq_struct, 0, "el_int", tmpmx);
    int16_t *eli = (int16_t *)mxGetData(tmpmx);
    tmpmx = mxCreateNumericMatrix(1, num_pulses, mxINT16_CLASS, mxREAL);
    mxSetField(iq_struct, 0, "az_int", tmpmx);
    int16_t *azi = (int16_t *)mxGetData(tmpmx);
    
	i = 0;  // index for MATLAB array
	j = 0;  // index for reader's packet count
	ip = 0;
	
	while (!feof(fid) && run && i<num_pulses) {
		fread(&ph, sizeof(struct packet_header), 1, fid);
		//mexPrintf("%c/%d  ", ph.type, ph.size);
		if (ph.size > (uint16_t)sizeof(PXIQPulse)) {
			mexPrintf("Probably a corrupted file. ph.size:%d > sizeof(PXIQPulse):%d\n",
					ph.size, (uint16_t)sizeof(PXIQPulse));
			break;
		}

		switch (ph.type) {
			case 'p':

				read_size = MIN(ph.size, sizeof(pulse));

				fread(&pulse, sizeof(PXIQPulseHeader), 1, fid);
				fread(&pulse.X[0][0].i, sizeof(PXIQ), pulse.ngate, fid);
				fread(&pulse.X[1][0].i, sizeof(PXIQ), pulse.ngate, fid);

				pulse_size = sizeof(PXIQPulseHeader) + 2*pulse.ngate*sizeof(PXIQ);

				if (j % 100 == 0 || pulse.vm & PXIQ_MARKER_SWEEP_BEGIN || pulse.vm & PXIQ_MARKER_SWEEP_END) {
				//if (1) {
					struct timeval time = {pulse.time_sec, pulse.time_usec};
					datestr(msg, &time, 32, 1);
					mexPrintf("%s %08d EL:% 5.2f (%+6.3f) AZ:%6.2f (%+6.2f) g%d v%d %05d/%05d\n",
						msg, pulse.n, pulse.el_deg, pulse.vel_dps, pulse.az_deg, pulse.vaz_dps, pulse.ngate, pulse.vm,
							j, i);
				}
				if (pulse.vm != PXIQ_MARKER_NULL || (nrhs > 2)) {

					float
						*c1i = &ch1_i[i*ngate],
						*c1q = &ch1_q[i*ngate],
						*c2i = &ch2_i[i*ngate],
						*c2q = &ch2_q[i*ngate];

                    int16_t *x1 = &pulse.X[0][0].i;
					int16_t *x2 = &pulse.X[1][0].i;

					for (k=0; k<pulse.ngate; k++) {
						*c1i++ = (float)*x1++;
						*c1q++ = (float)*x1++;
						*c2i++ = (float)*x2++;
						*c2q++ = (float)*x2++;
					}
					az[i] = pulse.az_deg;
					el[i] = pulse.el_deg;

                    eli[i] = pulse.el;
                    azi[i] = pulse.az;
                    
					i++;
				}
				//if (PXIQ_is_sweep_begin(pulse)) {
				//}

				break;
			default:
				mexPrintf("Ignored packet\n");
				break;
		}
		if ((ph.type == 'p' && i >= num_pulses) || (ph.type == 'p' && pulse.vm & PXIQ_MARKER_SWEEP_END))
			break;
		if (read_size < (int)ph.size) {
			mexPrintf("Seek forward %d B\n", ph.size-read_size);
			fseek(fid, ph.size-read_size, SEEK_CUR);
		}
		j++; // pulse in file counter;
	}
	fclose(fid);

	mxSetN(ch1, i);
	mxSetN(ch2, i);
	mxSetN(elmx, i);
	mxSetN(azmx, i);
    mxSetN(mxGetField(iq_struct, 0, "el_int"), i);
    mxSetN(mxGetField(iq_struct, 0, "az_int"), i);
}
