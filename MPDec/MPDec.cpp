// MPDec.cpp : Demo of message passing algorithm
// 
// Note
// rxLLR: ln(p(rx|0) / p(rx|1))
//
// Author: Alan Bao Jian ZHOU
// Email: zhoubj1990@gmail.com
// Created: 20140722
// Last-modified: 20140727


#include "stdafx.h"
#include "MPDecoder.h"


int _tmain(int argc, _TCHAR* argv[])
{
	// Start
	printf("----Message Passing Decoding Algorithm----\n");
	printf("Received LLR values: rxLLR[i] = ln(p(rx[i]|0) / p(rx[i]|1))\n");
	printf("Example:\n");
	printf("Channel code: (3, 1) binary repetition code\n");
	printf("Number of bits: nBit = 3\n");
	printf("Number of checks: nCheck = 2\n");
	printf("Number of edges: nEdge = 4\n");
	printf("Linear indices of edges: idxLinear = 0 1 2 5\n");
	printf("Modulation: BPSK (0 -> +1, 1 -> -1)\n");
	printf("PSD of AWGN: N0 = 1\n");
	printf("Received waveform: rx = 1 1 -1\n");
	printf("Received LLR values: rxLLR = 4 / N0 * rx = 4 4 -4\n");
	printf("Maximum number of iterations: nIterationMax = 10\n");
	printf("Estimated codeword: cHat = 0 0 0\n");
	printf("Number of iterations executed: nIteration = 1\n\n");

	// Code structure
	int nBit = 0;
	int nCheck = 0;
	int nEdge = 0;

	printf("nBit = ");
	scanf("%d", &nBit);

	printf("nCheck = ");
	scanf("%d", &nCheck);

	printf("nEdge = ");
	scanf("%d", &nEdge);

	double *idxLinear = new double[nEdge];

	printf("idxLinear = ");
	int iEdge = 0;
	for(iEdge = 0; iEdge < nEdge; iEdge++)
	{
		scanf("%lf", &idxLinear[iEdge]);
	}

	// MPDecoder
	MPDecoder *mpd = new MPDecoder(nBit, nCheck, nEdge, idxLinear);

	// Received LLR values, and decoding settings
	double *rxLLR = new double[nBit];
	int nIterationMax = 0;

	// Result
	double *cHat = new double[nBit];
	double *nIteration = new double;

	// Decoding loop
	// There's no input verification, so be cautious: the programm may crash due to improper inputs.
	while(1)
	{
		// Data
		printf("rxLLR = ");
		int iBit = 0;
		for(iBit = 0; iBit < nBit; iBit++)
		{
			scanf("%lf", &rxLLR[iBit]);
		}

		printf("nIterationMax = ");
		scanf("%d", &nIterationMax);

		// MP Decoding
		mpd->Decode(cHat, nIteration, rxLLR, nIterationMax);

		// Decoding result
		printf("cHat =");
		for(iBit = 0; iBit < nBit; iBit++)
		{
			printf(" %1.0lf", cHat[iBit]);
		}
		printf("\nnIteration = %1.0lf\n\n", *nIteration);

		// Exit or continue
		while(getchar() != '\n'); // Clear residue chars
		printf("Press ENTER to continue, or enter any other key to exit: ");
		if(getchar() != '\n')
		{
			break;
		}
	} // while Decoding loop

	// Release memory
	delete[] idxLinear;
	delete[] cHat;
	delete nIteration;
	delete mpd;
	delete[] rxLLR;

	// End
	return 0;
} // main