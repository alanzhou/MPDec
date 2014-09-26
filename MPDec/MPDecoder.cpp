// MPDecoder.cpp: Implementation of the MPDecoder class
// MPDecoder implements the message passing (MP) algorithm
// for decoding binary LDPC codes.
//
// Author: Alan Bao Jian ZHOU
// Email: zhoubj1990@gmail.com
// Created: 20140722
// Last-modified: 20140730

// Note
// i: traverse variable; idx: indexing variable.

#include "stdafx.h" // Only for MSVC++, remove otherwise.
#include <stdlib.h> // new, delete
#include <string.h> // memcpy, memset
#include <math.h> // tanh

#include "MPDecoder.h"

/*
MPDecoder::MPDecoder
Constructor, build structure and allocate memory.
Arguments
int nBit - input, number of bit nodes (number of columns of the parity check matrix)
int nCheck - input, number of check nodes (number of rows of the parity check matrix)
int nEdge -	input, number of edges (number of 1's in the parity check matrix)
double *idxLinear - input, linear indices of edges: idxLinear[iEdge] = nCheck * idxBit[iEdge] + idxCheck[iEdge]
	memory managed by caller
Caution
idxLinear[iEdge] counts from ZERO in C/C++ here, not from ONE as in MATLAB.
*/
MPDecoder::MPDecoder(int nBit, int nCheck, int nEdge, double *idxLinear)
{
	// Accept input
	this->nBit = nBit;
	this->nCheck = nCheck;
	this->nEdge = nEdge;

	// Initialization
	CalculateDegree(idxLinear); // Calculate node degrees
	BuildEdge(); // Build edges
	InitMessage(); // Initialize messages
} // MPDecoder


/*
MPDecoder::~MPDecoder
Destructor, release memory.
*/
MPDecoder::~MPDecoder()
{
	// Delete memory of nodes
	delete[] idxBit;
	delete[] idxCheck;
	delete[] degBit;
	delete[] degCheck;
	
	// Delete memory of edges
	for(int iBit = 0; iBit != nBit; ++iBit)
	{
		delete[] idxEdge2Bit[iBit];
	}
	delete[] idxEdge2Bit;

	for(int iCheck = 0; iCheck != nCheck; ++iCheck)
	{
		delete[] idxEdge2Check[iCheck];
	}
	delete[] idxEdge2Check;

	// Delete momory of messages
	delete[] msgChannel;
	delete[] msgBit;
	delete[] msgCheck;
	delete[] msgBit2Check;
	delete[] msgCheck2Bit;
} // ~MPDecoder


/*
MPDecoder::Decode
Perform MP decoding.
Arguments
double *cHat - output, estimated codeword, length: nBit, memory managed by caller
double *nIteration - output, number of iterations executed, memory managed by caller
double *rxLLR - input, pointer to the array of the received LLR values
int nIterationMax - input, maximum number of iterations
*/
void MPDecoder::Decode(double *cHat, double *nIteration, double *rxLLR, int nIterationMax)
{
	// Update channel messages, and set as initial bit messages
	UpdateChannel(rxLLR);

	HDD(cHat);// Hard decision
	if(IsValid(cHat))
	{
		*nIteration = 0; // Reset nIteration
		return;
	}

	// Important: Clear check-to-bit messages
	memset(msgCheck2Bit, 0, sizeof(double) * nEdge);

	// Decoding iterations
	for(*nIteration = 0; *nIteration != nIterationMax; ++(*nIteration))
	{
		UpdateCheck(); // Check node update
		UpdateBit(); // Bit node update

		HDD(cHat); // Hard decision
		if(IsValid(cHat))
		{
			++(*nIteration); // Increment nIteration
			return;
		}
	} // for *nIteration
} // Decode


/*
MPDecoder::CalculateDegree
Calculate degrees of bit nodes and check nodes.
Arguments
double *idxLinear - input, linear indices of edges: idxLinear[iEdge] = nCheck * idxBit[iEdge] + idxCheck[iEdge]
	memory managed by caller
*/
void MPDecoder::CalculateDegree(double *idxLinear)
{
	// Allocate memory for indices and degrees
	idxBit = new int[nEdge];
	idxCheck = new int[nEdge];
	degBit = new int[nBit];
	degCheck = new int[nCheck];

	memset(idxBit, -1, sizeof(int) * nEdge);
	memset(idxCheck, -1, sizeof(int) * nEdge);
	memset(degBit, 0, sizeof(int) * nBit);
	memset(degCheck, 0, sizeof(int) * nCheck);

	// Calculate indices and degrees
	// Caution: idxLinear[iEdge] = nCheck * idxBit[iEdge] + idxCheck[iEdge]
	// Note: Consider optimizing the code by introducing temporary variables like
	//			idxBit_iEdge = idxBit[iEdge];
	//		inside the loop to avoid frequent and redundant memory access.
	//		Similarly, do optimization for other loops.
	for(int iEdge = 0; iEdge != nEdge; ++iEdge)
	{
		// Caution: Since idxLinear_MAX = 64800^2 ~ 2^31.9674, when doing type
		// casting, use unsigned int, not int.
		idxBit[iEdge] = (unsigned int)idxLinear[iEdge] / nCheck;
		idxCheck[iEdge] = (unsigned int)idxLinear[iEdge] % nCheck;
		++degBit[ idxBit[iEdge] ];
		++degCheck[ idxCheck[iEdge] ];
	}
} // CalculateDegree


/*
MPDecoder::BuildEdge
Build edges.
*/
void MPDecoder::BuildEdge()
{
	// Allocate memory for edges
	idxEdge2Bit = new int *[nBit];
	memset(idxEdge2Bit, NULL, sizeof(int *) * nBit);
	for(int iBit = 0; iBit != nBit; ++iBit)
	{
		idxEdge2Bit[iBit] = new int[ degBit[iBit] ];
		memset(idxEdge2Bit[iBit], -1, sizeof(int) * degBit[iBit]);
	}

	idxEdge2Check = new int *[nCheck];
	memset(idxEdge2Check, NULL, sizeof(int *) * nCheck);
	for(int iCheck = 0; iCheck != nCheck; ++iCheck)
	{
		idxEdge2Check[iCheck] = new int[ degCheck[iCheck] ];
		memset(idxEdge2Check[iCheck], -1, sizeof(int) * degCheck[iCheck]);
	}

	// Build edges
	int *iEdge2Bit = new int[nBit];
	int *iEdge2Check = new int[nCheck];

	memset(iEdge2Bit, 0, sizeof(int) * nBit);
	memset(iEdge2Check, 0, sizeof(int) * nCheck);

	for(int iEdge = 0; iEdge != nEdge; ++iEdge)
	{
        int idxBit_iEdge = idxBit[iEdge];
		idxEdge2Bit[ idxBit_iEdge ][ iEdge2Bit[ idxBit_iEdge ] ] = iEdge;
		++iEdge2Bit[ idxBit_iEdge ];

        int idxCheck_iEdge = idxCheck[iEdge];
		idxEdge2Check[ idxCheck_iEdge ][ iEdge2Check[ idxCheck_iEdge ] ] = iEdge;
		++iEdge2Check[ idxCheck_iEdge ];
	}

	delete[] iEdge2Bit;
	delete[] iEdge2Check;
} // BuildEdge


/*
MPDecoder::InitMessage
Initialize messages.
*/
void MPDecoder::InitMessage()
{
	// Allocate memory for messages
	msgChannel = new double[nBit];
	msgBit = new double[nBit];
	msgCheck = new double[nCheck];
	msgBit2Check = new double[nEdge];
	msgCheck2Bit = new double[nEdge];

	memset(msgChannel, 0, sizeof(double) * nBit);
	memset(msgBit, 0, sizeof(double) * nBit);
	memset(msgCheck, 0, sizeof(double) * nCheck);
	memset(msgBit2Check, 0, sizeof(double) * nEdge);
	memset(msgCheck2Bit, 0, sizeof(double) * nEdge);
} // InitMessage


/*
MPDecoder::UpdateChannel
Update channel messages, and set as initial bit messages
Arguments
double *rxLLR - input, received LLR values
*/
void MPDecoder::UpdateChannel(double *rxLLR)
{
	// Update channel messages
	memcpy(msgChannel, rxLLR, sizeof(double) * nBit);
	// Set channel messages as initial bit messages
	memcpy(msgBit, msgChannel, sizeof(double) * nBit);
} // UpdateChannel


/*
MPDecoder::UpdateCheck
Update check nodes.
*/
void MPDecoder::UpdateCheck()
{
	// Update bit-to-check messages, and transform into tanh(msg/2) field
	for(int iEdge = 0; iEdge != nEdge; ++iEdge)
	{
		msgBit2Check[iEdge] = tanh((msgBit[ idxBit[iEdge] ] - msgCheck2Bit[iEdge]) / 2);
	}

	// Update check messages, in tanh(msg/2) field
	// Caution: memset works on a byte-by-byte basis only, and memset(p, 1, size_t) is not working!
	for(int iCheck = 0; iCheck != nCheck; ++iCheck)
	{
        double msgCheck_iCheck = 1;
        int degCheck_iCheck = degCheck[iCheck];
        int *idxEdge2Check_iCheck = idxEdge2Check[iCheck];
		for(int iEdge2Check = 0; iEdge2Check != degCheck_iCheck; ++iEdge2Check)
		{
			msgCheck_iCheck *= msgBit2Check[ idxEdge2Check_iCheck[iEdge2Check] ];
		} // for iEdge2Check
        msgCheck[iCheck] = msgCheck_iCheck;
	} // for iCheck
} // UpdateCheck


/*
MPDecoder::UpdateBit
Update bit nodes.
*/
void MPDecoder::UpdateBit()
{
	// Update check-to-bit messages, exit from tanh(msg/2) field
	for(int iEdge = 0; iEdge != nEdge; ++iEdge)
	{
		// Caution: Division by 0
		msgCheck2Bit[iEdge] = 2 * atanh(msgCheck[ idxCheck[iEdge] ] / msgBit2Check[iEdge]);
	}

	// Update bit messages
	for(int iBit = 0; iBit != nBit; ++iBit)
	{
        double msgBit_iBit = msgChannel[iBit];
        int degBit_iBit = degBit[iBit];
        int *idxEdge2Bit_iBit = idxEdge2Bit[iBit];
		for(int iEdge2Bit = 0; iEdge2Bit != degBit_iBit; ++iEdge2Bit)
		{
			msgBit_iBit += msgCheck2Bit[ idxEdge2Bit_iBit[iEdge2Bit] ];
		} // for iEdge2Bit
        msgBit[iBit] = msgBit_iBit;
	} // for iBit
} // UpdateBit


/*
MPDecoder::HDD
Set cHat with hard decision decoding (HDD).
Arguments
double *cHat - output, estimated codeword, length: nBit, memory managed by caller
*/
void MPDecoder::HDD(double *cHat)
{
    // Use memset to initialize cHat as all-0
    // memset(cHat, 0, sizeof(double) * nBit);
	for(int iBit = 0; iBit != nBit; ++iBit)
	{
		if(msgBit[iBit] < 0)
		{
			cHat[iBit] = 1;
		}
        else
        {
            cHat[iBit] = 0;
        }
	}
} // HDD


/*
MPDecoder::IsValid
Check if cHat is valid by examine the parity check condition.
Arguments
double *cHat - output, estimated codeword, length: nBit, memory managed by caller
Return
bool - true if cHat is valid; otherwise, false.
*/
bool MPDecoder::IsValid(double *cHat)
{
	for(int iCheck = 0; iCheck != nCheck; ++iCheck)
	{
		int parity = 0; // The parity check result for each check

        int degCheck_iCheck = degCheck[iCheck];
        int *idxEdge2Check_iCheck = idxEdge2Check[iCheck];
		for(int iEdge2Check = 0; iEdge2Check != degCheck_iCheck; ++iEdge2Check)
		{
			parity += (int)cHat[ idxBit[ idxEdge2Check_iCheck[iEdge2Check] ] ];
		} // for iBit

		if(parity & 1)
		{
			return false;
		}
	} // for iCheck

	return true;
} // IsValid