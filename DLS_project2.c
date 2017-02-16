/*
 * DLS_Project2.c
 *
 *  Created on: May 5, 2014
 *      Author: Paresh B
 */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>

#define PI     3.14159265358979323846
#define MAXDOF  10	/* maximum number of degrees of freedom for robot */
#define MAXDOFPLUS  MAXDOF+1 /* maximum number of robot coordinate frames */

typedef struct
{
        float x;
        float y;
} vector;


/* GLOBALS */

	/* file pointers */
	FILE *fopen(), *fp_armfile , *fp_trajfile , *fp_anglesfile;


	/* arm parameters */
	int nmbrdof;		/* number of degrees of freedom */
	int nmbrofpos;
	float LAMBDA2;
	float maxtheta;		/* maximum change in configuration between frames*/
	float lngth[MAXDOF];	/* link lengths */
	float theta[MAXDOF];	/* joint angle positions */
	float thetadot[MAXDOF]; /* joint angle velocities */
	//vector link_frame[MAXDOFPLUS]; /* world position (x,y) of individual */
				       /* link coordinate frames and end effector */

	vector efjacob[MAXDOF]; /* Jacobian for end effector */

	vector x_desired;	/* desired end effector trajectory position */
	vector xdot_desired;	/* desired end effector trajectory velocity */
	vector xdot_commanded;	/* commanded end effector velocity */

	vector arm[MAXDOF];	/* screen positions for displaying */
				/* the robot arm */
	float x_actual;
	float y_actual;
	float x_position[200];
	float y_position[200];
	float delx[100];
	float dely[100];
	float alpha[4];
	float delxi=200;
	float delyi=200;
	int temp_flag=0;
	float t1, t2;
	int flag1 = 0;
	float des_pos_x=0;
	float des_pos_y=0;

defnarm() /* read in description of robot arm */
{
	int i;
	fscanf(fp_armfile,"%d",&nmbrdof);
	if (nmbrdof > MAXDOF) printf("ERROR: TOO MANY DEGREES OF FREEDOM IN ROBOT\n");
	printf("%%\tNumber of Links: %d\n",nmbrdof);
	fscanf(fp_armfile,"%f", &maxtheta);
	printf("%%\tMaximum delta theta: %f\n", maxtheta);
	printf("%%\tLink Lengths \tAngles (degrees)\n");

	for (i = 0; i < nmbrdof; i++)
	{
		fscanf(fp_armfile,"%f%f",&lngth[i],&theta[i]);
		printf("%%\t%f\t%f\n",lngth[i],theta[i]);
		//theta[i]=theta[i]*PI/180; /* convert to radians */
	}
}


efjacobian() /* Calculate end effector Jacobian */
{
	int i,j;
	float psi; /* absolute joint angle */
	vector aj[MAXDOF]; /* absolute angle jacobian */
	/* calculate terms of Jacobian in absolute angles */
	psi=0;
	for (i = 0; i < nmbrdof; i++)
	{
		psi=psi+alpha[i];
		aj[i].x=-lngth[i]*sin(psi);
		aj[i].y= lngth[i]*cos(psi);
	}

	/* Jacobian in relative angles is a sum of absolute terms */
	for (i = 0; i < nmbrdof; i++)
	{
		efjacob[i].x=0;
		efjacob[i].y=0;
		for (j = i; j < nmbrdof; j++)
		{
			efjacob[i].x=efjacob[i].x+aj[j].x;
			efjacob[i].y=efjacob[i].y+aj[j].y;
			printf("Jacobian values %f %f\n",efjacob[i].x,efjacob[i].y);
		}
	}
	if (temp_flag == 1)
	{
		flag1++;
		solve();
	}
}

solve() /* solve the inverse kinematics using damped least squares */
{
	int i;
	float determ,temp;
	float t8=0;
	float t6=0;
	float t7=0;
	vector efpseudo[MAXDOF];
	float jjT[2][2];
	/* calculate the Jacobian times the Jacobian transpose (jjT) */
	jjT[0][0]=0;
	jjT[1][0]=0;
	jjT[0][1]=0;
	jjT[1][1]=0;
	for (i = 0; i < nmbrdof; i++)
	{
		jjT[0][0]=jjT[0][0]+efjacob[i].x*efjacob[i].x;
		jjT[1][0]=jjT[1][0]+efjacob[i].y*efjacob[i].x;

		jjT[0][1]=jjT[0][1]+efjacob[i].x*efjacob[i].y;
		jjT[1][1]=jjT[1][1]+efjacob[i].y*efjacob[i].y;
	}
	printf("Jacob JJt is %f %f %f %f\n",jjT[0][0],jjT[0][1],jjT[1][0],jjT[1][1]);
	/* add damping factor to the diagonal */
	jjT[0][0]=jjT[0][0]+LAMBDA2;
	jjT[1][1]=jjT[1][1]+LAMBDA2;
	printf("Jacob JJt + lambda sq I is %f %f %f %f\n",jjT[0][0],jjT[0][1],jjT[1][0],jjT[1][1]);

	/* calculate the inverse of jjT */
	determ=jjT[0][0]*jjT[1][1]-jjT[0][1]*jjT[1][0];
	temp=jjT[0][0];
	jjT[0][0]=jjT[1][1]/determ;
	jjT[0][1]=-jjT[0][1]/determ;
	jjT[1][0]=-jjT[1][0]/determ;
	jjT[1][1]=temp/determ;
	printf("Inverse of Jacob JJt is %f %f %f %f\n",jjT[0][0],jjT[0][1],jjT[1][0],jjT[1][1]);

	/* multiply by jT for damped least squares inverse, i.e. jT(jjT)^-1 */
	for (i = 0; i < nmbrdof; i++)
	{
		efpseudo[i].x=0;
		efpseudo[i].y=0;
		efpseudo[i].x=efpseudo[i].x+efjacob[i].x*jjT[0][0]
					   +efjacob[i].y*jjT[1][0];
		efpseudo[i].y=efpseudo[i].y+efjacob[i].x*jjT[0][1]
					   +efjacob[i].y*jjT[1][1];
	}
	for(i=0;i<nmbrdof;i++)
		printf("Pseudo inverse values are %f %f\n",efpseudo[i].x,efpseudo[i].y);

	/* multiply Jacobian inverse times the commanded end effector */
	/* velocity to obtain the desired joint angle velocity */
	for (i = 0; i < nmbrdof; i++)
	{
	thetadot[i]=efpseudo[i].x*t1+efpseudo[i].y*t2;
	}
	/* integrate the joint velocity to find joint position */
	i=0;
	while(i < nmbrdof)
	{
		alpha[i]=alpha[i]+thetadot[i];
		i++;
	}
	float sum_angles = 0;
	float square_angles=0;
	x_actual = 0;
	y_actual=0;
	i=0;
	while(i < nmbrdof)
	{
		sum_angles = sum_angles + alpha[i];
		square_angles = square_angles + alpha[i]*alpha[i];
		x_actual = x_actual + 10*cos(sum_angles);
		y_actual = y_actual + 10*sin(sum_angles);
		i++;
	}

	t8=delxi;
	t6=delyi;
	t7=sqrt(t8*t8 + t6*t6);
	delxi = des_pos_x - x_actual;
	delyi = des_pos_y - y_actual;
	float temp4 = 0;
	float temp5 = 0;
	temp4 =  (sqrt((delxi*delxi) + (delyi*delyi)));
	temp5 = sqrt(square_angles);
	if ((temp4<t7) && (temp5<maxtheta))
	{
	t1 = delxi;
	t2 = delyi;
	temp_flag = 1;
	efjacobian();
	}
	else
	{
		temp_flag = 0;
		x_actual = 0;
		y_actual=0;
		sum_angles=0;
		square_angles=0;
	}

}

calculate_diff()
{
	int i;
	float velocity_norm;
	float velocity_norm_temp=100;
	for (i=0; i<(nmbrofpos-1); i++)
	{
		delx[i] = x_position[i+1]-x_position[i];
		dely[i] = y_position[i+1]-y_position[i];
	}
	for(i=0 ; i<(nmbrofpos-1); i++)
	{
		t1 = delx[i];
		t2 = dely[i];
		velocity_norm=sqrt(t1*t1 + t2*t2);
		if(velocity_norm < velocity_norm_temp)
			velocity_norm_temp=velocity_norm;
	}
		float t1=0;
		float t2=0;
		t1 = maxtheta*2;
		t2=velocity_norm_temp/t1;
		LAMBDA2 = t2*t2;
}

int main()
{
	int i;
fp_armfile = fopen("arm","r"); 					//open the file to read
fp_trajfile = fopen("trajectory","r");	//open the file to read
defnarm(); /* read in arm description */
i=0;
while(i<nmbrdof)
{
		alpha[i]=theta[i];
		i++;

}
		//defntraj();
int j;
	fscanf(fp_trajfile,"%d",&nmbrofpos);
	printf("%%\tNumber of desired positions: %d\n",nmbrofpos);
	printf("%%\tDesired x position \Desired y position\n");
	for (j = 0; j < (nmbrofpos); j++)
	{
		fscanf(fp_trajfile,"%f%f",&x_position[j],&y_position[j]);
		printf("%%\t%f\t%f\n",x_position[j],y_position[j]);
	}
calculate_diff();
int k;
int m=0;
fp_anglesfile = fopen("angles","w");		//open file for writing
for(i=0;i<nmbrdof;i++)
	fprintf(fp_anglesfile,"%f ",theta[i]);
fprintf(fp_anglesfile,"\n");
for(k=0;k<nmbrofpos-1;k++)
{
	m++;
	efjacobian(); /* calculate end effector Jacobian */
t1=delx[k];
t2=dely[k];
des_pos_x=x_position[m];
des_pos_y=y_position[m];
solve(); /* solve for new joint angles */
delxi=100;
delyi=100;
	//open file for writing
flag1 = 0;
int z;
for(z=0;z<nmbrdof;z++)
	fprintf(fp_anglesfile,"%f ",alpha[z]);
fprintf(fp_anglesfile,"\n");
}
}


