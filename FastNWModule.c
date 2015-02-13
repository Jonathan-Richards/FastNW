/**********************************************************************
* Author: Jonathan Richards
* Date: 1/1/2015
*
* Module containing a stripped-down implementation of Needleman–Wunsch
* algorithm that only returns the best global alignment score, not any
* alignments. Allows for gap extension penalites. Space scaling lowered
* to O(min(n, m)) using the observation that only the current and
* previous row of the score matrix must be stored.
**********************************************************************/

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

__inline int mymax(int a, int b) {
  return a > b ? a : b;
}

int withExtend(const char *shorter, const char *longer, size_t width, size_t height,
	int match, int mismatch, int gap, int gap_extend) {
	//current and previous row given not in a gap, in a down gap, and in a right gap
	int *cur;
	int *prev;
	int *cur_down;
	int *prev_down;
	int *cur_right;
	int *prev_right;

	int *temp; //for switching cur and prev

	//loop variables
	size_t i;
	size_t j;

	//return value
	int ret;

	//allocate rows of matrix
	cur = malloc(width*sizeof(int));
	prev = malloc(width*sizeof(int));
	cur_right = malloc(width*sizeof(int));
	prev_right = malloc(width*sizeof(int));
	cur_down = malloc(width*sizeof(int));
	prev_down = malloc(width*sizeof(int));

	/*************** Initial assignment of cur ***************/
	cur[0] = 0;
	cur_right[0] = INT_MIN/2;
	cur_down[0] = INT_MIN/2;
	for (i=1; i<width; i++) {
		cur[i] = INT_MIN/2;
		cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
		cur_down[i] = INT_MIN/2;
	}

	/***************** Assign rest of matrix *****************/
	for (j=1; j<height; j++) {
		//current becomes previous
		temp = prev;
		prev = cur;
		cur = temp;

		temp = prev_down;
		prev_down = cur_down;
		cur_down = temp;

		temp = prev_right;
		prev_right = cur_right;
		cur_right = temp;

		cur[0] = INT_MIN/2;
		cur_right[0] = INT_MIN/2;
		cur_down[0] = mymax(prev[0]+gap, prev_down[0]+gap_extend);

		//calculate current row
		for (i=1; i<width; i++) {
			
			//calculate score after diagonal path
			if (shorter[i-1] == longer[j-1]) {
				cur[i] = mymax(prev[i-1], mymax(prev_right[i-1], prev_down[i-1])) + match;
			} else {
				cur[i] = mymax(prev[i-1], mymax(prev_right[i-1], prev_down[i-1])) + mismatch;
			}

			//calculate score after downward path
			cur_down[i] = mymax(prev[i] + gap, prev_down[i] + gap_extend);

			//calculate score after rightward path
			cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
		}
	}

	ret = mymax(cur[width-1], mymax(cur_right[width-1], cur_down[width-1]));
	free(prev);
	free(prev_down);
	free(cur);
	free(cur_right);
	free(cur_down);
	return ret;
}

static PyObject * NWScore(PyObject *self, PyObject *args) {
	//input strings
	const char *shorter;
	const char *longer;

	const char *temp; //for switching shorter and longer

	//penalty scores
	int match;
	int mismatch;
	int gap;
	int gap_extend = INT_MIN;

	//dimensions of matrix
	size_t width;
	size_t height;

	int ret; //return value

	if (!PyArg_ParseTuple(args, "ssiii|i", &shorter, &longer, &match, &mismatch, &gap, &gap_extend))
		return NULL;

	//find shorter and longer inputs
	if (strlen(shorter) > strlen(longer)) {
		temp = longer;
		longer = shorter;
		shorter = temp;
	}

	//assign dimensions of matrix
	width = strlen(shorter)+1;
	height = strlen(longer)+1;

	//call appropriate function
	if (gap_extend == INT_MIN) {
		gap_extend = gap;
	}
	ret = withExtend(shorter, longer, width, height,
		match, mismatch, gap, gap_extend);

	return Py_BuildValue("i", ret);
}

static PyMethodDef NWMethods[] = {
    {"score",  NWScore, METH_VARARGS,
     "Compute a Needleman–Wunsch score"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initFastNW(void) {
    (void) Py_InitModule("FastNW", NWMethods);
}

int main(int argc, char *argv[]) {
	Py_SetProgramName(argv[0]);
	Py_Initialize();
	initFastNW();
	return 0;
}

/*
int main(int argc, char *argv[]) {
	//input strings
	char *shorter;
	char *longer;

	//penalty scores
	int match = 1;
	int mismatch = -5;
	int gap = -10;
	int gap_extend = -1;

	//dimensions of matrix
	size_t width;
	size_t height;

	//overwrite default penalty scores
	if (argc > 3) {
		match = atoi(argv[3]);
		mismatch = atoi(argv[4]);
		gap = atoi(argv[5]);
	}
	if (argc > 6) {
		gap_extend = atoi(argv[6]);
	}
	
	//find shorter and longer inputs
	if (strlen(argv[1]) < strlen(argv[2])) {
		shorter = argv[1];
		longer = argv[2];
	} else {
		longer = argv[1];
		shorter = argv[2];
	}

	//assign dimensions of matrix
	width = strlen(shorter)+1;
	height = strlen(longer)+1;

	if (argc > 6 || argc <= 3) {
		printf("%d, %d, %d\n", width, height, withExtend(shorter, longer, width, height,
			match, mismatch, gap, gap_extend));
	} else {
		printf("%d\n", noExtend(shorter, longer, width, height,
			match, mismatch, gap));
	}

	return 0;
}
*/
