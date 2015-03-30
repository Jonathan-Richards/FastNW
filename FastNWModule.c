/**********************************************************************
* FastNW: Fast Needleman-Wunsch
* Copyright (C) 2014 Jonathan Richards
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
* 
* Written by Jonathan Richards, jonrds@gmail.com
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

typedef enum {false, true} bool;

typedef struct {
	char *shorter;
	char *longer;
	int match;
	int mismatch;
	int gap;
	int gap_extend;
	bool switched;
} Arguments;

const Arguments FAILED = {
	NULL, NULL, 0, 0, 0, 0, false
};

typedef struct {
	int score;
	char *align1;
	char *align2;
} Alignment;

typedef struct {
	int score;
	size_t index;
} HirschReturn;

const HirschReturn NEED_MEM = {
	0, -1
};

typedef struct {
	int *cur;
	int *cur_right;
	int *cur_down;
} ScoreReturn;

const ScoreReturn NO_MEM = {
	0, NULL, NULL
};

typedef enum {NONE, DOWN, RIGHT, ANY} Direction;

typedef struct {
	int index;
	Direction left;
	Direction right;
} PartitionReturn;

//for quickly counting a Needleman Wunsch _score_
//returns a ScoreReturn object that must be freed after use
//returning entire bottom row allows use in Partition
ScoreReturn Score(const char *horizontal, size_t hl, size_t hr,
	const char *vertical, size_t vl, size_t vr,
	int match, int mismatch, int gap, int gap_extend,
	Direction start_direction) {

	//dimensions of matrix
	size_t width = hr-hl+1;
	size_t height = vr-vl+1;

	//loop variables
	size_t i;
	size_t j;

	//return value
	ScoreReturn ret;

	//current and previous row given not in a gap, in a down gap, and in a right gap
	int *cur = malloc(width*sizeof(int));
	int *prev = malloc(width*sizeof(int));
	int *cur_right = malloc(width*sizeof(int));
	int *prev_right = malloc(width*sizeof(int));
	int *cur_down = malloc(width*sizeof(int));
	int *prev_down = malloc(width*sizeof(int));
	int *temp; //for switching cur and prev


	/******************** Check Memory ***********************/
	if (cur==NULL || prev==NULL || cur_right==NULL
		|| prev_right==NULL || cur_down==NULL || prev_down==NULL) {

		free(cur);
		free(prev);
		free(cur_right);
		free(prev_right);
		free(cur_down);
		free(prev_down);
		return NO_MEM;
	}

	/*************** Initial assignment of cur ***************/
	cur[0] = 0;
	cur_right[0] = INT_MIN/4;
	cur_down[0] = INT_MIN/4;
	for (i=1; i<width; i++) {
		cur[i] = INT_MIN/4;
		cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
		cur_down[i] = INT_MIN/4;
	}

	/******** Second row depends on start_direction **********/
	if (height > 1) {
		temp = prev;
		prev = cur;
		cur = temp;

		temp = prev_down;
		prev_down = cur_down;
		cur_down = temp;

		temp = prev_right;
		prev_right = cur_right;
		cur_right = temp;

		cur[0] = INT_MIN/4;
		cur_right[0] = INT_MIN/4;
		switch (start_direction) {
			case NONE : //cant use prev_right or cur_down
				cur_down[0] = INT_MIN/4;
				for (i=1; i<width; i++) {
					if (horizontal[hl+i-1] == vertical[vl])
						cur[i] = prev[i-1]+match;
					else
						cur[i] = prev[i-1]+mismatch;
					cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
					cur_down[i] = INT_MIN/4;
				}
				break;
			case DOWN : //can only use cur_down
				cur_down[0] = gap;
				for (i=1; i<width; i++) {
					cur[i] = INT_MIN/4;
					cur_right[i] = INT_MIN/4;
					cur_down[i] = INT_MIN/4;
				}
				break;
			case RIGHT : //cant use prev or cur_down
				cur_down[0] = INT_MIN/4;
				for (i=1; i<width; i++) {
					if (horizontal[hl+i-1] == vertical[vl])
						cur[i] = prev_right[i-1]+match;
					else
						cur[i] = prev_right[i-1]+mismatch;
					cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
					cur_down[i] = INT_MIN/4;
				}
				break;
			case ANY : //can use arrays as normal
				cur_down[0] = gap;
				for (i=1; i<width; i++) {
					if (horizontal[hl+i-1] == vertical[vl])
						cur[i] = mymax(prev[i-1], prev_right[i-1])+match;
					else
						cur[i] = mymax(prev[i-1], prev_right[i-1])+mismatch;
					cur_right[i] = mymax(cur[i-1] + gap, cur_right[i-1] + gap_extend);
					cur_down[i] = INT_MIN/4;
				}
				break;
			default :
				printf("ERROR: Initial direction error\n");
				break;
		}
	}

	/***************** Assign rest of matrix *****************/
	for (j=2; j<height; j++) {
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

		cur[0] = INT_MIN/4;
		cur_right[0] = INT_MIN/4;
		cur_down[0] = mymax(prev[0]+gap, prev_down[0]+gap_extend);

		//calculate current row
		for (i=1; i<width; i++) {
			
			//calculate score after diagonal path
			if (horizontal[hl+i-1] == vertical[vl+j-1]) {
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

	free(prev);
	free(prev_right);
	free(prev_down);

	ret.cur = cur;
	ret.cur_right = cur_right;
	ret.cur_down = cur_down;

	return ret;
}

//Full Needleman Wunsch algorithm, with added capability
//for starting and ending requirements (allowing it to be used with Hirsch)
HirschReturn NeedlemanWunsch(char *Z, char *W, size_t Z_spot,
	const char *horizontal, size_t hl, size_t hr,
	const char *vertical, size_t vl, size_t vr,
	int match, int mismatch, int gap, int gap_extend,
	Direction start_direction, Direction end_direction) {

	//for indexing
	size_t i;
	size_t j;

	//matrix dimensions
	size_t width = hr-hl+1;
	size_t height = vr-vl+1;

	//temporary calculations
	int from;
	int from_right;
	int from_down;

	//score matrix
	int *mat = malloc(width*height*sizeof(int));
	int *mat_right = malloc(width*height*sizeof(int));
	int *mat_down = malloc(width*height*sizeof(int));

	//0=none, 1=right, 2=down
	int trace;
	int *mat_dir = malloc(width*height*sizeof(int));
	int *mat_right_dir = malloc(width*height*sizeof(int));
	int *mat_down_dir = malloc(width*height*sizeof(int));

	//backwards alignments
	size_t rev_spot;
	char *rev_Z = malloc((width+height)*sizeof(char));
	char *rev_W = malloc((width+height)*sizeof(char));

	HirschReturn ret;

	if (mat==NULL || mat_right==NULL || mat_down==NULL
		|| mat_dir==NULL || mat_right_dir==NULL || mat_down_dir==NULL
		|| rev_Z==NULL || rev_W==NULL) {

		free(mat);
		free(mat_right);
		free(mat_down);
		free(mat_dir);
		free(mat_right_dir);
		free(mat_down_dir);
		free(rev_Z);
		free(rev_W);
		return NEED_MEM;
	}

/*
	printf(horizontal);
	printf("\n");
	printf(vertical);
	printf("\n");

	printf("width = %d\n", width);
	printf("height = %d\n", height);
*/

	/********************** First row ************************/
	mat[0] = 0;
	mat_dir[0] = -1;

	mat_right[0] = INT_MIN/4;
	mat_right_dir[0] = -1;
	
	mat_down[0] = INT_MIN/4;
	mat_down_dir[0] = -1;
	
	for (i=1; i<width; i++) {
		mat[i] = INT_MIN/4;
		mat_dir[i] = -1;

		from = mat[i-1] + gap;
		from_right = mat_right[i-1] + gap_extend;
		if (from > from_right) {
			mat_right[i] = from;
			mat_right_dir[i] = 0;
		} else {
			mat_right[i] = from_right;
			mat_right_dir[i] = 1;
		}

		mat_down[i] = INT_MIN/4;
		mat_down_dir[i] = -1; //added down
	}

	/******** Second row depends on start_direction **********/
	if (height > 1) {
		j = width;

		mat[j] = INT_MIN/4;
		mat_dir[j] = -1;

		mat_right[j] = INT_MIN/4;
		mat_right_dir[j] = -1;

		switch (start_direction) {
			case NONE : //cant use prev_right or cur_down
				mat_down[j] = INT_MIN/4;
				mat_down_dir[j] = -1;

				for (i=1; i<width; i++) {
					if (horizontal[hl+i-1] == vertical[vl])
						mat[j+i] = mat[i-1]+match;
					else
						mat[j+i] = mat[i-1]+mismatch;
					mat_dir[j+i] = 0;

					from = mat[j+i-1] + gap;
					from_right = mat_right[j+i-1] + gap_extend;
					if (from > from_right) {
						mat_right[j+i] = from;
						mat_right_dir[j+i] = 0;
					} else {
						mat_right[j+i] = from_right;
						mat_right_dir[j+i] = 1;
					}
					
					mat_down[j+i] = INT_MIN/4;
					mat_down_dir[j+i] = -1;
				}
				break;
			case DOWN : //can only use cur_down
				//printf("Starting down\n");
				mat_down[j] = gap;
				mat_down_dir[j] = 0;

				for (i=1; i<width; i++) {
					mat[j+i] = INT_MIN/4;
					mat_dir[j+i] = -1;

					mat_right[j+i] = INT_MIN/4;
					mat_right_dir[j+i] = -1;

					mat_down[j+i] = INT_MIN/4;
					mat_down_dir[j+i] = -1;
				}
				break;
			case RIGHT : //cant use prev or cur_down
				mat_down[j] = INT_MIN/4;
				mat_down_dir[j] = -1;

				for (i=1; i<width; i++) {
					if (horizontal[hl+i-1] == vertical[vl])
						mat[j+i] = mat_right[i-1]+match;
					else
						mat[j+i] = mat_right[i-1]+mismatch;
					mat_dir[j+i] = 0;

					from = mat[j+i-1] + gap;
					from_right = mat_right[j+i-1] + gap_extend;
					if (from > from_right) {
						mat_right[j+i] = from;
						mat_right_dir[j+i] = 0;
					} else {
						mat_right[j+i] = from_right;
						mat_right_dir[j+i] = 1;
					}
					
					mat_down[j+i] = INT_MIN/4;
					mat_down_dir[j+i] = -1;
				}
				break;
			case ANY : //can use arrays as normal
				mat_down[j] = gap;
				mat_down_dir[j] = 0;

				for (i=1; i<width; i++) {
					from = mat[i-1];
					from_right = mat_right[i-1];
					if (from > from_right) {
						mat[j+i] = from;
						mat_dir[j+i] = 0;
					} else {
						mat[j+i] = from_right;
						mat_dir[j+i] = 1;
					}
					if (horizontal[hl+i-1] == vertical[vl])
						mat[j+i] += match;
					else
						mat[j+i] += mismatch;

					from = mat[j+i-1] + gap; //added j
					from_right = mat_right[j+i-1] + gap_extend; //added j
					if (from > from_right) {
						mat_right[j+i] = from;
						mat_right_dir[j+i] = 0;
					} else {
						mat_right[j+i] = from_right;
						mat_right_dir[j+i] = 1;
					}

					mat_down[j+i] = INT_MIN/4;
					mat_down_dir[j+i] = -1;
				}
				break;
			default :
				printf("ERROR: Initial direction error\n");
				break;
		}
	}

	/***************** Assign rest of matrix *****************/
	for (j=2*width; j<width*height; j+=width) {

		mat[j] = INT_MIN/4;
		mat_dir[j] = -1;

		mat_right[j] = INT_MIN/4;
		mat_right_dir[j] = -1;

		from = mat[j-width]+gap;
		from_down = mat_down[j-width]+gap_extend;
		if (from > from_down) {
			mat_down[j] = from;
			mat_down_dir[j] = 0;
		} else {
			mat_down[j] = from_down;
			mat_down_dir[j] = 2;
		}

		//calculate current row
		for (i=1; i<width; i++) {
			
			//calculate score after diagonal path
			from = mat[j-width+i-1];
			from_right = mat_right[j-width+i-1];
			from_down = mat_down[j-width+i-1];
			if (from > from_right && from > from_down) {
				mat[j+i] = from;
				mat_dir[j+i] = 0;
			} else if (from_right > from_down) {
				mat[j+i] = from_right;
				mat_dir[j+i] = 1;
			} else {
				mat[j+i] = from_down;
				mat_dir[j+i] = 2;
			}
			if (horizontal[hl+i-1] == vertical[vl+j/width-1]) {
				mat[j+i] += match;
			} else {
				mat[j+i] += mismatch;
			}

			//calculate score after rightward path
			from = mat[j+i-1] + gap; //removed -width
			from_right = mat_right[j+i-1] + gap_extend; //removed -width
			if (from > from_right) {
				mat_right[j+i] = from;
				mat_right_dir[j+i] = 0;
			} else {
				mat_right[j+i] = from_right;
				mat_right_dir[j+i] = 1;
			}

			//calculate score after downward path
			from = mat[j-width+i] + gap;
			from_down = mat_down[j-width+i] + gap_extend;
			if (from > from_down) {
				mat_down[j+i] = from;
				mat_down_dir[j+i] = 0;
			} else {
				mat_down[j+i] = from_down;
				mat_down_dir[j+i] = 2;
			}
		}
	}

	/*
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat[j+i]);
		}
		printf("\n");
	}
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat_dir[j+i]);
		}
		printf("\n");
	}
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat_right[j+i]);
		}
		printf("\n");
	}
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat_right_dir[j+i]);
		}
		printf("\n");
	}
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat_down[j+i]);
		}
		printf("\n");
	}
	for (j=0; j<width*height; j+=width) {
		for (i=0; i<width; i++) {
			printf("%d, ", mat_down_dir[j+i]);
		}
		printf("\n");
	}
	*/
	

	/*********** Matrix completed, begin backtrace **************/

	//calculate backtrace starting position
	i = width*height-1;
	switch (end_direction) {
		case NONE :
			ret.score = mat[i];
			trace = 0;
			//rev_Z[0] = horizontal[hr-1];
			//rev_W[0] = vertical[vr-1];
			break;
		case RIGHT :
			ret.score = mat_right[i];
			trace = 1;
			//rev_Z[0] = horizontal[hr-1];
			//rev_W[0] = '-';
			break;
		case DOWN :
			//printf("start trace down\n");
			ret.score = mat_down[i];
			trace = 2;
			//rev_Z[0] = '-';
			//rev_W[0] = vertical[vr-1];
			break;
		case ANY :
			if (mat[i] > mat_right[i] && mat[i] > mat_down[i]) {
				ret.score = mat[i];
				trace = 0;
				//rev_Z[0] = horizontal[hr-1];
				//rev_W[0] = vertical[vr-1];
			} else if (mat_right[i] > mat_down[i]) {
				ret.score = mat_right[i];
				trace = 1;
				//rev_Z[0] = horizontal[hr-1];
				//rev_W[0] = '-';
			} else {
				//printf("hi\n");
				ret.score = mat_down[i];
				trace = 2;
				//rev_Z[0] = '-';
				//rev_W[0] = vertical[vr-1];
				//printf("Trace = %d\n", trace);
			}
			break;
		default :
			printf("ERROR: end direction error\n");
			break;
	}

	//printf("Populating\n");

	//populate rev alignments
	j = height-1;
	i = width-1;
	rev_spot = 0;
	while (j > 0 || i > 0) {
		//printf("i=%d, j=%d, t=%d\n", i, j, trace);
		switch (trace) {
			case 0 :
				trace = mat_dir[j*width+i];
				i--;
				j--;
				rev_Z[rev_spot] = horizontal[hl+i];
				rev_W[rev_spot] = vertical[vl+j];
				break;
			case 1 :
				trace = mat_right_dir[j*width+i];
				i--;
				rev_Z[rev_spot] = horizontal[hl+i];
				rev_W[rev_spot] = '-';
				break;
			case 2 :
				trace = mat_down_dir[j*width+i];
				j--;
				rev_Z[rev_spot] = '-';
				rev_W[rev_spot] = vertical[vl+j];
				break;
			default :
				printf("ERROR: traceback error\n");
				exit(0);
				break;
		}
		
		rev_spot++;
	}

	ret.index = Z_spot + rev_spot;

	/************ Put rev alignments into Z and W ************/

	rev_Z[rev_spot] = '\0';
	rev_W[rev_spot] = '\0';

/*
	printf("Aligning ");
	printf(rev_Z);
	printf(" and ");
	printf(rev_W);
	printf("\n");
*/

	for (rev_spot; rev_spot > 0; rev_spot--, Z_spot++) {
		//printf("%d, %d\n", rev_spot, Z_spot);
		Z[Z_spot] = rev_Z[rev_spot-1];
		W[Z_spot] = rev_W[rev_spot-1];
	}
	if (Z_spot != ret.index) {
		printf("ERROR: indexing problem\n");
	}

	free(mat);
	free(mat_right);
	free(mat_down);
	free(mat_dir);
	free(mat_right_dir);
	free(mat_down_dir);
	free(rev_Z);
	free(rev_W);

	//printf("Done\n");

	return ret;

}

PartitionReturn Partition(ScoreReturn ScoreL, ScoreReturn ScoreR, size_t width,
	int gap, int gap_extend) {
	size_t i;
	size_t j=width;
	int best = INT_MIN;
	int score;
	PartitionReturn ret;

	/*
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreL.cur[i]);
	}
	j=width;
	printf("\n");
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreL.cur_down[i]);
	}
	j=width;
	printf("\n");
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreL.cur_right[i]);
	}
	j=width;
	printf("\n");
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreR.cur[i]);
	}
	j=width;
	printf("\n");
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreR.cur_down[i]);
	}
	j=width;
	printf("\n");
	for (i=0; i<=width; i++,j--) {
		printf("%d, ", ScoreR.cur_right[i]);
	}
	j=width;
	printf("\n");
	*/

	//need some way to force other partitions into ending in gap or not
	for (i=0; i<=width; i++, j--) {
		//neither gap
		score = ScoreL.cur[i] + ScoreR.cur[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = NONE;
			ret.right = NONE;
		}

		//both down gaps. Have to correct scores
		score = ScoreL.cur_down[i] + ScoreR.cur_down[j] - gap + gap_extend;
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = DOWN;
			ret.right = DOWN;
		}

		//one is down, other is not
		score = ScoreL.cur[i] + ScoreR.cur_down[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = NONE;
			ret.right = DOWN;
		}

		//one is down, other is not
		score = ScoreL.cur_down[i] + ScoreR.cur[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = DOWN;
			ret.right = NONE;
		}

		//any combination of right and cur can be represented here
		score = ScoreL.cur_right[i] + ScoreR.cur[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = RIGHT;
			ret.right = NONE;
		}

		/*
		score = ScoreL.cur_right[i] + ScoreR.cur_right[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = RIGHT;
			ret.right = RIGHT;
		}

		score = ScoreL.cur[i] + ScoreR.cur_right[j];
		if (score > best) {
			best = score;
			ret.index = i;
			ret.left = NONE;
			ret.right = RIGHT;
		}
		*/
	}
	//printf("score = %d, %d\n", best, ret.index);

	free(ScoreL.cur);
	free(ScoreL.cur_right);
	free(ScoreL.cur_down);
	free(ScoreR.cur);
	free(ScoreR.cur_right);
	free(ScoreR.cur_down);
	
	return ret;
}

//recursive function for hirshberg algorithm
HirschReturn Hirsch(char *Z, char *W, size_t Z_spot,
	const char *horizontal, const char *rev_hor, size_t hl, size_t hr,
	const char *vertical, const char *rev_vert, size_t vl, size_t vr,
	int match, int mismatch, int gap, int gap_extend,
	Direction start_direction, Direction end_direction) {

	//get input string lengths
	size_t width = hr-hl;
	size_t height = vr-vl;

	//for indexing
	//size_t i;
	//size_t j;
	size_t h_mid;
	size_t v_mid;

	//alignment scores fo partitioning
	ScoreReturn ScoreL;
	ScoreReturn ScoreR;

	PartitionReturn pres;

	HirschReturn ret; //return value
	HirschReturn res; //result from NeedlemanWunsch

	ret.score = 0; //the relative score of this recursion call
	ret.index = Z_spot; //the absolute position in aligned strings
	//printf("hor: %d, %d\n", hl, hr);
	//printf("vert: %d, %d\n", vl, vr);
	//printf("width: %d\n", width);
	//printf("height: %d\n", height);

	if (width*height <= 1000000 || width==1 || height==1) {
		//printf("Args: %d, %d, %d, %d, %d\n", hl, hr, vl, vr, Z_spot);

		/*		
		printf("start = %d\n", start_direction);
		printf("end = %d\n", end_direction);

		for (i=hl; i<hr; i++)  {
			printf("%c", horizontal[i]);
		}
		printf("\n");
		for (j=vl; j<vr; j++)  {
			printf("%c", vertical[j]);
		}
		printf("\n");
		*/

		res = NeedlemanWunsch(Z, W, Z_spot, horizontal, hl, hr,
			vertical, vl, vr,
			match, mismatch, gap, gap_extend,
			start_direction, end_direction);

		ret.score = res.score;
		ret.index = res.index;
	} else {
		v_mid = (vl+vr)/2; //split vertical in half

		ScoreL = Score(horizontal, hl, hr,
			vertical, vl, v_mid,
			match, mismatch, gap, gap_extend, start_direction);
		ScoreR = Score(rev_hor, strlen(rev_hor)-hr, strlen(rev_hor)-hl,
			rev_vert, strlen(rev_vert)-vr, strlen(rev_vert)-v_mid,
			match, mismatch, gap, gap_extend, end_direction);

		//partition horizontal
		pres = Partition(ScoreL, ScoreR, width, gap, gap_extend);
		h_mid = hl+pres.index;
		//printf("pres.left = %d\n", pres.left);
		//printf("pres.right = %d\n", pres.right);
		
		/*
		for (i=hl; i<hr; i++)  {
			printf("%c", horizontal[i]);
		}
		printf("\n");
		printf("%d\n", h_mid);

		for (j=vl; j<vr; j++)  {
			printf("%c", vertical[j]);
		}
		printf("\n");
		*/
		/*
		printf("hor: %d, %d, %d\n", hl, h_mid, hr);
		printf("vert: %d, %d, %d\n", vl, v_mid, vr);
		printf(rev_hor);
		printf("\n");
		printf(rev_vert);
		printf("\n");
		*/

		res = Hirsch(Z, W, Z_spot,
			horizontal, rev_hor, hl, h_mid,
			vertical, rev_vert, vl, v_mid,
			match, mismatch, gap, gap_extend,
			start_direction, pres.left);
		ret.score = res.score;
		Z_spot = res.index;

		res = Hirsch(Z, W, Z_spot,
			horizontal, rev_hor, h_mid, hr,
			vertical, rev_vert, v_mid, vr,
			match, mismatch, gap, gap_extend,
			pres.right, end_direction);
		ret.score += res.score;
		ret.index = res.index;

		//have to remember that this is actually 1 gap, not 2
		if (pres.left == DOWN && pres.right == DOWN)
			ret.score += gap_extend-gap;
	}

	return ret;

}

//interprets python arguments
Arguments GetArguments(PyObject *args) {
	char *temp; //for switching longer and shorter
	Arguments arguments; //return value
	arguments.gap_extend = INT_MIN;
	
	//parse python args
	if (!PyArg_ParseTuple(args, "ssiii|i",
		&arguments.shorter, &arguments.longer, &arguments.match,
		&arguments.mismatch, &arguments.gap, &arguments.gap_extend))
		return FAILED;

	//find shorter and longer inputs
	if (arguments.switched = (strlen(arguments.shorter) > strlen(arguments.longer))) {
		temp = arguments.shorter;
		arguments.shorter = arguments.longer;
		arguments.longer = temp;
	}

	//if gap_extend isn't specified, must be equal to gap
	if (arguments.gap_extend == INT_MIN) {
		arguments.gap_extend = arguments.gap;
	}

	return arguments;
}

//handler for score method from python
static PyObject * NWScore(PyObject *self, PyObject *args) {
	ScoreReturn res; //return value
	int width;
	int ret;

	Arguments arguments = GetArguments(args);
	if (!arguments.shorter)
		return NULL;
	width = strlen(arguments.shorter);

	res = Score(arguments.shorter, 0, strlen(arguments.shorter),
		arguments.longer, 0, strlen(arguments.longer),
		arguments.match, arguments.mismatch, arguments.gap, arguments.gap_extend,
		ANY);

	ret = mymax(res.cur[width], mymax(res.cur_right[width], res.cur_down[width]));
	free(res.cur);
	free(res.cur_right);
	free(res.cur_down);

	return Py_BuildValue("i", ret);
}

//handler for align method from python
static PyObject * Align(PyObject *self, PyObject *args) {
	HirschReturn res; //result from hirschberg algorithm
	PyObject *ret; //return value

	char *Z; //short alignment
	char *W; //long alignment

	char *rev_hor; //reverse of input strings
	char *rev_vert;

	//input string sizes
	size_t width;
	size_t height;

	size_t i;

	Arguments arguments = GetArguments(args);
	if (!arguments.shorter)
		return NULL;

	width = strlen(arguments.shorter);
	height = strlen(arguments.longer);
	rev_hor = malloc((width+1)*sizeof(char));
	rev_vert = malloc((height+1)*sizeof(char)); // extra +1 for \0 to use strlen
	Z = malloc((width+height)*sizeof(char));
	W = malloc((width+height)*sizeof(char));

	for (i=0; i<width; i++) {
		rev_hor[width-i-1] = arguments.shorter[i];
	}
	rev_hor[width] = '\0';
	for (i=0; i<height; i++) {
		rev_vert[height-i-1] = arguments.longer[i];
	}
	rev_vert[height] = '\0';


	res = Hirsch(Z, W, 0,
		arguments.shorter, rev_hor, 0, width,
		arguments.longer, rev_vert, 0, height,
		arguments.match, arguments.mismatch, arguments.gap, arguments.gap_extend,
		ANY, ANY);

	Z[res.index] = '\0';
	W[res.index] = '\0';

	if (arguments.switched)
		ret = Py_BuildValue("[s,s,i]", W, Z, res.score);
	else
		ret = Py_BuildValue("[s,s,i]", Z, W, res.score);

	free(Z);
	free(W);

	return ret;
}

//handler for qalign method from python
static PyObject * QAlign(PyObject *self, PyObject *args) {
	HirschReturn res; //result from hirschberg algorithm
	PyObject *ret; //return value

	char *Z; //short alignment
	char *W; //long alignment

	//input string sizes
	size_t width;
	size_t height;

	Arguments arguments = GetArguments(args);
	if (!arguments.shorter)
		return NULL;

	width = strlen(arguments.shorter);
	height = strlen(arguments.longer);
	Z = malloc((width+height)*sizeof(char));
	W = malloc((width+height)*sizeof(char));

	res = NeedlemanWunsch(Z, W, 0,
		arguments.shorter, 0, width,
		arguments.longer, 0, height,
		arguments.match, arguments.mismatch, arguments.gap, arguments.gap_extend,
		ANY, ANY);

	//printf("Done2\n");

	Z[res.index] = '\0';
	W[res.index] = '\0';

	if (arguments.switched)
		ret = Py_BuildValue("[s,s,i]", W, Z, res.score);
	else
		ret = Py_BuildValue("[s,s,i]", Z, W, res.score);

	free(Z);
	free(W);

	return ret;
}

static PyMethodDef NWMethods[] = {
    {"score",  NWScore, METH_VARARGS,
     "Compute a Needleman–Wunsch score"},
    {"align", Align, METH_VARARGS,
	 "Compute a Needleman-Wunsch alignment using the Hirschberg Algorithm"},
	{"qalign", QAlign, METH_VARARGS,
	 "Force a Needleman-Wunsch alignment"},
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
