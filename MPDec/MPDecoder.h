// MPDecoder.h: Definition of the MPDecoder class
// MPDecoder implements the message passing (MP) algorithm
// for decoding binary LDPC codes.
//
// Author: Alan Bao Jian ZHOU
// Email: zhoubj1990@gmail.com
// Created: 20140722
// Last-modified: 20140727


#ifndef _MP_DECODER_H
#define _MP_DECODER_H

#include <math.h> // log

class MPDecoder
{
	// Caution: Interfacing to MATLAB MEX requires variables being double!
	/* Constructor and destructor */
public:
	MPDecoder(int nBit, int nCheck, int nEdge, double *idxLinear); // Build structure and allocate memory
	~MPDecoder(); // Release memory

	/* Member variables */
	// Code structure
private:
	int nBit; // Number of bit nodes (number of columns of the parity check matrix)
	int nCheck; // Number of check nodes (number of rows of the parity check matrix)
	int nEdge; // Number of edges (number of 1's in the parity check matrix)

private:
	int *idxBit; // Bit indices of edges
	int *idxCheck; // Check indices of edges
	int *degBit; // Bit node degrees
	int *degCheck; // Check node degrees
	int **idxEdge2Bit; // idxEdge2Bit[iBit][ degBit[iBit] ]: edge indices of edges to bit node iBit
	int **idxEdge2Check; // idxEdge2Check[iCheck][ degCheck[iCheck] ]: edge indices of edges to check node iCheck

	// Messages
private:
	double *msgChannel; // Channel messages
	double *msgBit; // Bit messages
	double *msgCheck; // Check messages
	double *msgBit2Check; // Bit-to-check messages
	double *msgCheck2Bit; // Check-to-bit messages

	/* Member functions */
	// Interface
public:
	void Decode(double *cHat, double *nIteration, double *rxLLR, int nIterationMax); // Perform MP decoding

	// Building-block functions
private:
	void CalculateDegree(double *idxLinear); // Calculate degrees of bit nodes and check nodes
	void BuildEdge(); // Build edges
	void InitMessage(); // Initialize messages

	void UpdateChannel(double *rxLLR); // Update channel messages
	void UpdateCheck(); // Update check nodes
	void UpdateBit(); // Update bit nodes

	void HDD(double *cHat); // Hard decision decoding (HDD)
	bool IsValid(double *cHat); // Check if cHat is valid
}; // class MPDecoder

// Before MSVC++ 11, atanh is not included in MS version of math.h
inline double atanh(double x)
{
    // atanh(1) and atanh(-1) are set to be 19.07 and -19.07
    // to avoid infinite numbers being used. For details, check
    // http://www.mathworks.com/help/comm/ref/ldpcdecoder.html
    if(x <= -1)
    {
        return -19.07;
    }
    else if (x >= 1)
    {
        return 19.07;
    }
    else
    {
        return 0.5 * log((1+x) / (1-x));
    }
} // atanh

#endif	// _MP_DECODER_H