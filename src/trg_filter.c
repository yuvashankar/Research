 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

#include <stdlib.h>
#include <stdio.h>



int num_odd,num_regular;

int num_res_reg_right,num_res_reg_wrong, num_res_odd_right,num_res_odd_wrong;


int mode; //mode is 1 for regular, 2 for odd 

int last_presented; // (1 for regular, 2 for odd);
int responded; // We only allow one response to a stimulus


int num_since_odd;


float time;
long val1;
long val2;

float last_time;


float av_reg_res_time_right;
float av_reg_res_time_wrong;
float av_odd_res_time_right;
float av_odd_res_time_wrong;



int read_trigger(FILE* in_file)
{
	int ret_code;
	ret_code=fscanf(in_file,"%f",&time);
	ret_code=fscanf(in_file,"%ld",&val1);
	ret_code=fscanf(in_file,"%ld",&val2);
	return(ret_code);
}

int process_file(FILE* in_file, FILE* out_file)
{
	int ret_code=1;
	char buff[64];

	int code,button;

	// Rad the first 3 lines, they are nto important for us 
	fgets(buff,64,in_file);
	fgets(buff,64,in_file);
	fgets(buff,64,in_file);


	//Now line by line 3 numbers in each
	while(ret_code==1){
		ret_code=read_trigger(in_file);
		// Check, there could be a --- in the last collunm
		if(ret_code==0){
			printf("ERROR \n");
			exit(0);
		}else{
			code = val2 & 255;
			button=(val2 >>8)&3;

			
			if((code==1) && (mode==0)){
				num_regular++;
				mode=1;
				last_presented=1;
				last_time=time;
				responded=0;
				if(num_since_odd!=-1) num_since_odd++;
				fprintf(out_file,"\n %d\t REGULAR\t %d \t ",
					num_odd+num_regular, num_since_odd);
			}
			if((code==2) && (mode==0)){
				num_odd++;
				mode=2;
				last_presented=2;
				last_time=time;
				responded=0;
				num_since_odd=0;
				fprintf(out_file,"\n %d \t ODD    \t %d  \t ",
					num_odd+num_regular, num_since_odd);
			}
			if(code==0){
				mode=0;
			}

			// NOW WE CHECK THE BUTTONS
			if((button==1) && (responded==0)){
				responded=1;
				if(last_presented==1){
					num_res_reg_right++;
					float r_t=time-last_time;
					av_reg_res_time_right+=r_t;
					fprintf(out_file," RIGHT \t %f ",r_t);
				}
				if(last_presented==2){
					num_res_odd_wrong++;
					float r_t=time-last_time;
					av_odd_res_time_wrong+=r_t;
					fprintf(out_file," WRONG \t %f ",r_t);
				}
				
			}
			if((button==2) && (responded==0)){
				responded=1;
				if(last_presented==1){
					num_res_reg_wrong++;
					float r_t=time-last_time;
					av_reg_res_time_wrong+=r_t;
					fprintf(out_file," WRONG \t %f ",r_t);
				}
				if(last_presented==2){
					num_res_odd_right++;
					float r_t=time-last_time;
					av_odd_res_time_right+=r_t;
					fprintf(out_file," RIGHT \t %f",r_t);
				}
			}
			//fprintf(out_file," %d ",val2);
			//fprintf(out_file," Condition  %d Button %d \n",code,button);

			

		}
		

		
	}
	printf("There were %d experiments, %d regular %d odd \n",
		num_regular+num_odd,num_regular,num_odd);

	printf("Responses \n");
	printf("Regular  %d right %d wrong \n",num_res_reg_right,num_res_reg_wrong);
	printf("Odd      %d right %d wrong \n",num_res_odd_right,num_res_odd_wrong);

	printf("Response Times \n");
	printf("Regular average time right %f wrong %f \n",
			av_reg_res_time_right/((float)(num_res_reg_right)),
			av_reg_res_time_wrong/((float)(num_res_reg_wrong)));
	printf("Odd average time: right %f wrong %f \n",
			av_odd_res_time_right/((float)(num_res_odd_right)),
			av_odd_res_time_wrong/((float)(num_res_odd_wrong)));
	

	return(0);
	

}



int main(int nargs, char* args[])
{
	FILE * filter_file;
	FILE* in_file;
	FILE* out_file;


	if(nargs!=3){
		printf("Usage: %s input_file.trg out_file \n",args[0]);
		return(1);
	}


	in_file=fopen(args[1],"r");
        if(in_file==NULL){
                printf("ERROR: input_file %s does not exist \n",args[1]);
                return(1);
        }

	out_file=fopen(args[2],"w");
	        if(out_file==NULL){
                printf("ERROR: input_file %s does not exist \n",args[2]);
                return(1);
        }


	//Initialize Gloabls
	num_odd=0;
	num_regular=0;
	num_res_reg_right=0;
	num_res_reg_wrong=0;
	num_res_odd_right=0;
	num_res_odd_wrong=0;
	last_presented=0;
	responded=0;

	
 	num_since_odd=-1;

	last_time=0.;
	av_reg_res_time_right=0.;
	av_reg_res_time_wrong=0.;
	av_odd_res_time_right=0.;
	av_odd_res_time_wrong=0.;


	printf(" \n Start Processing file:  %s  \n",args[1]);
	process_file(in_file,out_file);

	fclose(in_file);
	fclose(out_file);

	return(0);
}
