/*This code simulates system of hard circular disks using event driven method. Particles do elastic collisions with each other. Initial positions and velocities are chosen randomly. --------written by Dipanjan Mandal, dipkar.308@gmail.com*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mt19937ar.h"

#define L 50 //system size
#define N 60 // Tptal number of particles
#define T_MAX 10000//maximum event
#define max_v 20.0 // Maximum velocity of a particle
#define sigma 2.0 //particle radius
#define pair_dim (N*N/2-N/2)//total pair-collision events need to be calculated

double tm;
double pair_collision[pair_dim],wall_collision[N];
int first_in_pair[pair_dim],second_in_pair[pair_dim];
int wall_number[N];

typedef struct vector
{
	double x,y; 
}	vecr;

vecr pos[N],vel[N];

FILE *pipe;

void take_input()
{
	/*reads in seed for random number generator */
	long seedval;
	//printf("\nEnter value of seed : ");
	//scanf("%ld",&seedval);
	seedval=time(NULL);
	//seedval=989568;
	init_genrand(seedval);
}
int get_pair_index(int i, int j)
{
	int ind;

	if(j<i)
	{
		ind=i+N*j-(j+1)*(j+2)/2;
	}			

	if(j>i)
	{
		ind=j+N*i-(i+1)*(i+2)/2;
	}
	return ind;
}
void store_index_to_pair()
{
	int i,j,ind;
	
	ind=0;
	for(j=0;j<N;j++)
	{
		for(i=0;i<N;i++)
		{
			if(j<i)
			{
				first_in_pair[ind]=i;
				second_in_pair[ind]=j;
				ind++;
			}
		}
	}
}
double get_distance(int ind, double rx, double ry)
{
	double dist;
	dist=sqrt((pos[ind].x-rx)*(pos[ind].x-rx)+(pos[ind].y-ry)*(pos[ind].y-ry));
	return dist;
}
void initialize_pos_vel_randomly()
{
	int i,n,FLAG;
	double rx,ry,dist,d1,d2,d3,d4;

	n=0;
	while(n<N)
	{
		FLAG=0;
		rx=L*genrand_real3();
		ry=L*genrand_real3();	
		//check overlap with other particles
		if(n>0)
		{
			for(i=0;i<n;i++)
			{
				dist=get_distance(i,rx,ry);
				if(dist>=2.0*sigma)
				{
					FLAG=1;
				}
				else
				{
					FLAG=0;
					break;
				}
			}
		}
		else
		{
			FLAG=1;
		}
		//check overlap with the four walls
		if(FLAG==1)
		{
			d1=fabs(L-rx);
			d2=fabs(L-ry);
			d3=fabs(rx);
			d4=fabs(ry);
			if(d1>=sigma && d2>=sigma && d3>=sigma && d4>=sigma)
			{
				pos[n].x=rx;
				pos[n].y=ry;
				n++;
			}
		}
	}
	for(i=0;i<N;i++)
	{
		vel[i].x=2.0*max_v*genrand_real3();
		vel[i].y=2.0*max_v*genrand_real3();
		vel[i].x=vel[i].x-max_v;
		vel[i].y=vel[i].y-max_v;
	}
}
//calculate next collision time between pair particles k and l
double get_next_collision_time(int k, int l)
{
	double dist_x,dist_y,dv_x,dv_y;
	double delxdotdelv,delxdotdelx,delvdotdelv,gamma,t_pair;

	dist_x=pos[k].x-pos[l].x;
	dist_y=pos[k].y-pos[l].y;
	dv_x=vel[k].x-vel[l].x;
	dv_y=vel[k].y-vel[l].y;

	delxdotdelv=dist_x*dv_x+dist_y*dv_y;
	delxdotdelx=dist_x*dist_x+dist_y*dist_y;
	delvdotdelv=dv_x*dv_x+dv_y*dv_y;

	gamma=delxdotdelv*delxdotdelv-delvdotdelv*(delxdotdelx-4.0*sigma*sigma);
	if(gamma>0.0 && delxdotdelv<0.0)
	{
		t_pair=-(delxdotdelv+sqrt(gamma))/(1.0*delvdotdelv);
	}
	else
	{
		t_pair=1e6;//infinity
	}
	return t_pair;
}
//pair collision velocity update
void pair_collision_velocity_update(int k, int l, double velk[2], double vell[2])
{
	double dist_x,dist_y,dv_x,dv_y,unit_perp_x,unit_perp_y;

	dist_x=pos[k].x-pos[l].x;
	dist_y=pos[k].y-pos[l].y;
	dv_x=vel[k].x-vel[l].x;
	dv_y=vel[k].y-vel[l].y;

	unit_perp_x=dist_x/sqrt(dist_x*dist_x+dist_y*dist_y);
	unit_perp_y=dist_y/sqrt(dist_x*dist_x+dist_y*dist_y);

	velk[0]=vel[k].x-(dv_x*unit_perp_x+dv_y*unit_perp_y)*unit_perp_x;
	velk[1]=vel[k].y-(dv_x*unit_perp_x+dv_y*unit_perp_y)*unit_perp_y;

	vell[0]=vel[l].x+(dv_x*unit_perp_x+dv_y*unit_perp_y)*unit_perp_x;
	vell[1]=vel[l].y+(dv_x*unit_perp_x+dv_y*unit_perp_y)*unit_perp_y;
}
//wall collision velocity update
void wall_collision_velocity_update(int wall_index, int wall_no, double velk[2])
{
	if(wall_no==0)
	{
		velk[0]=-vel[wall_index].x;
		velk[1]=vel[wall_index].y;
	}
	else if(wall_no==1)
	{
		velk[0]=vel[wall_index].x;
		velk[1]=-vel[wall_index].y;
	}
	else if(wall_no==2)
	{
		velk[0]=-vel[wall_index].x;
		velk[1]=vel[wall_index].y;
	}
	else if(wall_no==3)
	{
		velk[0]=vel[wall_index].x;
		velk[1]=-vel[wall_index].y;
	}
}
//wall collision
double wall_collision_time(int k)
{
	double t,tp;

	if(vel[k].x>0.0)
	{
		if(vel[k].y>0.0)
		{
			t=(-pos[k].x+L-sigma)/(1.0*vel[k].x);
			wall_number[k]=0;
			tp=(-pos[k].y+L-sigma)/(1.0*vel[k].y);
			if(t>tp)
			{
				t=tp;
				wall_number[k]=1;
			}
		}
		else if(vel[k].y<0.0)
		{
			t=(-pos[k].x+L-sigma)/(1.0*vel[k].x);
			wall_number[k]=0;
			tp=(-pos[k].y+sigma)/(1.0*vel[k].y);
			if(t>tp)
			{
				t=tp;
				wall_number[k]=3;
			}
		}
		else if(vel[k].y==0.0)
		{
			t=(-pos[k].x+L-sigma)/(1.0*vel[k].x);
			wall_number[k]=0;
		}
	}
	else if(vel[k].x<0.0)
	{
		if(vel[k].y>0.0)
		{
			t=(-pos[k].x+sigma)/(1.0*vel[k].x);
			wall_number[k]=2;
			tp=(-pos[k].y+L-sigma)/(1.0*vel[k].y);
			if(t>tp)
			{
				t=tp;
				wall_number[k]=1;
			}
		}
		else if(vel[k].y<0.0)
		{
			t=(-pos[k].x+sigma)/(1.0*vel[k].x);
			wall_number[k]=2;
			tp=(-pos[k].y+sigma)/(1.0*vel[k].y);
			if(t>tp)
			{
				t=tp;
				wall_number[k]=3;
			}
		}
		else if(vel[k].y==0.0)
		{
			t=(-pos[k].x+sigma)/(1.0*vel[k].x);
			wall_number[k]=2;
		}
	}
	else if(vel[k].x==0.0)
	{
		if(vel[k].y>0.0)
		{
			t=(-pos[k].y+L-sigma)/(1.0*vel[k].y);
			wall_number[k]=1;
		}
		else if(vel[k].y<0.0)
		{
			t=(-pos[k].y+sigma)/(1.0*vel[k].y);
			wall_number[k]=3;
		}
		else if(vel[k].y==0.0)
		{
			t=1e6;
			printf("Velocity zero!!\n");
		}
	}
	return t;
}
//sort array in ascending order
int sort_array(int dim, double arr[])
{
	int i,j;
	double temp;
	int arr2[dim],temp2;

	for(i=0;i<dim;i++)
	{
		arr2[i]=i;
	}

	for(i=0;i<dim;i++)
	{
		for(j=i+1;j<dim;j++)
		{
			if(arr[i]>arr[j])
			{
				temp=arr[i];
				arr[i]=arr[j];
				arr[j]=temp;

				temp2=arr2[i];
				arr2[i]=arr2[j];
				arr2[j]=temp2;
			}
		}
	}
	return arr2[0];
}
void update_all_positions(double dt)
{
	int i;
	for(i=0;i<N;i++)
	{
		pos[i].x=pos[i].x+dt*vel[i].x;
		pos[i].y=pos[i].y+dt*vel[i].y;
	}
}
//evolve one collision event
void evolve_one_collision()
{
	int i,j,ind,ind_pair,ind_wall,k,l;
	double pair_time,wall_time,velk[2],vell[2];

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(j<i)
			{
				ind=get_pair_index(i,j);
				pair_collision[ind]=get_next_collision_time(i,j);
			}
		}
	}
	//get minimum pair-collision time
	ind_pair=sort_array(pair_dim,pair_collision);
	pair_time=pair_collision[0];
	k=first_in_pair[ind_pair];
	l=second_in_pair[ind_pair];

	for(i=0;i<N;i++)
	{
		wall_collision[i]=wall_collision_time(i);
	}
	//get minimum wall-collision time
	ind_wall=sort_array(N,wall_collision);
	wall_time=wall_collision[0];

	if(pair_time>wall_time)
	{
		update_all_positions(wall_time);
		wall_collision_velocity_update(ind_wall,wall_number[ind_wall],velk);
		vel[ind_wall].x=velk[0];
		vel[ind_wall].y=velk[1];
		tm+=wall_time;
	}
	else
	{
		update_all_positions(pair_time);
		pair_collision_velocity_update(k,l,velk,vell);
		vel[k].x=velk[0];
		vel[k].y=velk[1];
		vel[l].x=vell[0];
		vel[l].y=vell[1];
		tm+=pair_time;
	}
}
void call_gnuplot(int t)
{
	fprintf(pipe, "set title 'time=%.8e'\n",tm);
	fprintf(pipe,"p 'position.dat' u 2:3:4 w circle lw 2 lc rgb 'red' tit ''\n");
	fflush(pipe);
}
double calculate_total_energy(int time)
{
	int i;
	double s=0.0;
	for(i=0;i<N;i++)
	{
		s+=0.5*(vel[i].x*vel[i].x+vel[i].y*vel[i].y);
	}
	s=s/(1.0*N);
	return s;
}
int main()
{
	int i,j;
	char outfile[100];

	take_input();
	store_index_to_pair();
	initialize_pos_vel_randomly();

	pipe=popen("\\gnuplot","w");
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	fprintf(pipe, "set size square\n");
	fprintf(pipe, "unset key\n");
	fprintf(pipe, "set xrange[0:%d]\n",L);
	fprintf(pipe, "set yrange[0:%d]\n",L);

	tm=0.0;

	//evolve_one_collision();
	FILE *fp1;
	
	for(i=0;i<T_MAX;i++)
	{
		evolve_one_collision();
		if(i%20==0)
		{
			sprintf(outfile,"position.dat");
			fp1=fopen(outfile,"w");
  	 		fprintf(fp1,"#no\tx\ty\tsigma\n");
   		fclose(fp1);
			fp1=fopen(outfile,"a");	
			for(j=0;j<N;j++)
			{
				fprintf(fp1,"%d\t%.8e\t%.8e\t%.8e\n",j,pos[j].x,pos[j].y,sigma);
			}
			fclose(fp1);
			call_gnuplot(i);
			//printf("%.8e\t%.8e\n",tm,calculate_total_energy(i));
		}
	}
	return (0);
}
