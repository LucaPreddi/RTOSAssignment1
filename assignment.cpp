// Compile with: 
// $ g++ assignment.c -pthread -o assignment
// ./assignment

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>

//code of periodic tasks
void task1_code( );
void task2_code( );
void task3_code( );
void task4_code( );

double U;
//code of aperiodic tasks (if any)

//characteristic function of the thread, only for timing and synchronization
//periodic tasks

void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);

// Declaring global variables.

int T1T2 = 0, T1T4 = 0, T2T3 = 0;

// Declaring mutexes.

pthread_mutex_t T1T2_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t T1T4_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t T2T3_mutex = PTHREAD_MUTEX_INITIALIZER;

// Declaring time structs for tasks.

struct timespec start_1, end_1, start_2, end_2, start_3, end_3, start_4, end_4, start_5, end_5;
double Z12, Z14, blocking1, Z12a, Z23, blocking2, blocking3, blocking4;

// initialization of mutexes and conditions (only for aperiodic scheduling)

#define INNERLOOP 1000
#define OUTERLOOP 1000

#define NPERIODICTASKS 4
#define NAPERIODICTASKS 0
#define NTASKS NPERIODICTASKS + NAPERIODICTASKS

long int periods[NTASKS];
struct timespec next_arrival_time[NTASKS];
double WCET[NTASKS];
pthread_attr_t attributes[NTASKS];
pthread_t thread_id[NTASKS];
struct sched_param parameters[NTASKS];
int missed_deadlines[NTASKS];

void wastetime(){
	float guardie;
	for (int j = 0; j < OUTERLOOP; j++){
			guardie = rand()*rand();
	}
}

int main()
{
  	// set task periods in nanoseconds
	//the first task has period 200 millisecond
	//the second task has period 160 millisecond
	//the third task has period 100 millisecond
	//the third task has period 80 millisecond
	//you can already order them according to their priority; 
	//if not, you will need to sort them

  	periods[3] = 200000000;	// In nanoseconds.
  	periods[2] = 160000000;	// In nanoseconds.
  	periods[1] = 100000000;	// In nanoseconds.
  	periods[0] = 80000000;	// In nanoseconds.

  	//for aperiodic tasks we set the period equals to 0

	//this is not strictly necessary, but it is convenient to
	//assign a name to the maximum and the minimum priotity in the
	//system. We call them priomin and priomax.

  	struct sched_param priomax;
  	priomax.sched_priority = sched_get_priority_max(SCHED_FIFO);
  	struct sched_param priomin;
  	priomin.sched_priority = sched_get_priority_min(SCHED_FIFO);

	// set the maximum priority to the current thread (you are required to be
  	// superuser). Check that the main thread is executed with superuser privileges
	// before doing anything else.

  	if (getuid() == 0) pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);

  	// execute all tasks in standalone modality in order to measure execution times
  	// (use gettimeofday). Use the computed values to update the worst case execution
  	// time of each task.

 	int i;
  	for (i =0; i < NTASKS; i++){

		// initializa time_1 and time_2 required to read the clock

		struct timespec time_1, time_2;
		clock_gettime(CLOCK_REALTIME, &time_1);

		//we should execute each task more than one for computing the WCET
		//periodic tasks

 	    if (i==0) task1_code();
        if (i==1) task2_code();
      	if (i==2) task3_code();
		if (i==3) task4_code();
      		
      	//aperiodic tasks

		clock_gettime(CLOCK_REALTIME, &time_2);

		// compute the Worst Case Execution Time (in a real case, we should repeat this many times under
		//different conditions, in order to have reliable values

      		WCET[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec)
			       +(time_2.tv_nsec-time_1.tv_nsec);
      		printf("\nWorst Case Execution Time %d=%f \n", i, WCET[i]);
    	}

    // compute U
	

    // compute Ulub by considering the fact that we have harmonic relationships between periods
	double Ulub;
    	
	//if there are no harmonic relationships, use the following formula instead
	for(int i = 1 ; i<=NPERIODICTASKS; i++){
	
		Ulub = i*(pow(2.0,(1.0/i)) -1);

		switch(i){

			case(1):
				U = WCET[i-1]/periods[i-1] + blocking1/(periods[i-1]);	
				printf("U1: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(2):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + blocking2/(periods[i-1]);
				printf("U2: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(3):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + WCET[i-3]/periods[i-3] + blocking3/(periods[i-1]);
				printf("U3: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(4):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + WCET[i-3]/periods[i-3] + WCET[i-4]/periods[i-4] + blocking4/(periods[i-4]);
				printf("U4: %f, Ulub: %f\n", U, Ulub); fflush(stdout);

			break;
		}

		if (U > Ulub){
      		printf("\n U=%lf Ulub=%lf Non schedulable Task Set, %d", U, Ulub, i);
      		return(-1);
    	}
	}
	printf("\n U=%lf Ulub=%lf schedulable test", U, Ulub);
	fflush(stdout);
  	sleep(2);
	//check the sufficient conditions: if they are not satisfied, exit  


  	// set the minimum priority to the current thread: this is now required because 
	//we will assign higher priorities to periodic threads to be soon created
	//pthread_setschedparam

  	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  
  	// set the attributes of each task, including scheduling policy and priority
  	for (i = 0; i < NPERIODICTASKS; i++){

		//initializa the attribute structure of task i
      		pthread_attr_init(&(attributes[i]));

		//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		//not inherited from the main thread (pthread_attr_setinheritsched) 
      		pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

		//properly set the parameters to assign the priority inversely proportional 
		//to the period
      		parameters[i].sched_priority = priomin.sched_priority+NTASKS - i;

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
      		pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    }

    pthread_mutexattr_t mymutexattr; pthread_mutexattr_init(&mymutexattr); 
	pthread_mutexattr_setprotocol(&mymutexattr, PTHREAD_PRIO_PROTECT); 
	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority); 
	pthread_mutex_init(&T1T2_mutex, &mymutexattr); 
	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[1].sched_priority); 
	pthread_mutex_init(&T2T3_mutex, &mymutexattr);
	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority);
	pthread_mutex_init(&T1T4_mutex, &mymutexattr);

 	// aperiodic tasks

	//delare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	//declare variables to read the current time
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1);

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NPERIODICTASKS; i++)
    	{
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i];
		//then we compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}

	// create all threads(pthread_create)
  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
  	iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// join all threads (pthread_join)
  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
  	pthread_join( thread_id[3], NULL);

  	// Printing all the variables.
  	printf("\nThe variables are:\n");
  	printf("T1T2 = %d\n", T1T2);
  	printf("T1T4 = %d\n", T1T4);
  	printf("T2T3 = %d\n", T2T3);

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NTASKS; i++){
      	printf ("\nMissed Deadlines Task %d=%d", i, missed_deadlines[i]);
		fflush(stdout);
    }
  	exit(0);
}

// application specific task_1 code
void task1_code()
{
	//print the id of the current task
  	printf(" 1[ "); fflush(stdout);
	pthread_mutex_lock(&T1T2_mutex);
	clock_gettime(CLOCK_REALTIME, &start_1);
	wastetime();
	T1T2 += 1;
	clock_gettime(CLOCK_REALTIME, &end_1);
	pthread_mutex_unlock(&T1T2_mutex);

	pthread_mutex_lock(&T1T4_mutex);
	clock_gettime(CLOCK_REALTIME, &start_2);
	wastetime();
	T1T4 += 1;
	clock_gettime(CLOCK_REALTIME, &start_2);
	pthread_mutex_unlock(&T1T4_mutex);

	Z12 = 1000000000*(start_1.tv_sec - end_1.tv_sec) + (start_1.tv_nsec-end_1.tv_nsec);
	Z14 = 1000000000*(start_2.tv_sec - end_2.tv_sec) + (start_2.tv_nsec-end_2.tv_nsec);

	if (Z12>Z14){
		blocking1 = Z12; 
	}
	else if (Z14>Z12){
		blocking2 = Z14;
	}

  	//print the id of the current task
  	printf(".");
  	printf(" ]1 "); fflush(stdout);
}

//thread code for task_1 (used only for temporization)
void *task1( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

   	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i = 0;
  	for (i = 0; i < 100; i++)
    	{
      	// execute application specific code
		task1_code();

		// it would be nice to check if we missed a deadline here... why don't
		// you try by yourself?

		// sleep until the end of the current period (which is also the start of the
		// new one
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

		// the thread is ready and can compute the end of the current period for
		// the next iteration
 		
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task2_code()
{
	//print the id of the current task
  	printf(" 2[ "); fflush(stdout);

	pthread_mutex_lock(&T2T3_mutex);
	clock_gettime(CLOCK_REALTIME, &start_2);
	wastetime();
	T2T3 += 1;
	clock_gettime(CLOCK_REALTIME, &end_2);
	pthread_mutex_unlock(&T2T3_mutex);

	pthread_mutex_lock(&T1T2_mutex);
	clock_gettime(CLOCK_REALTIME, &start_3);
	wastetime();
	printf("%d", T1T2); fflush(stdout);
	clock_gettime(CLOCK_REALTIME, &end_3);
	pthread_mutex_unlock(&T1T2_mutex);

	//print the id of the current task
  	printf(" ]2 "); fflush(stdout);

  	Z12a = 1000000000*(start_2.tv_sec - end_2.tv_sec) + (start_2.tv_nsec-end_2.tv_nsec);
  	Z23 = 1000000000*(start_3.tv_sec - end_3.tv_sec) + (start_3.tv_nsec-end_3.tv_nsec);

  	if (Z12a>Z23){
  		blocking2 = Z12a;
  	}
  	else if(Z23>Z12a){
  		blocking2 = Z23;
  	}
}


void *task2( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i = 0;
  	for (i = 0; i < 100; i++)
    	{
      		task2_code();

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task3_code()
{
	//print the id of the current task
  	printf(" 3[ "); fflush(stdout);
  	pthread_mutex_lock(&T2T3_mutex);
  	clock_gettime(CLOCK_REALTIME, &start_4);
  	wastetime();
	printf("%d", T2T3); 
	clock_gettime(CLOCK_REALTIME, &end_4);
	pthread_mutex_unlock(&T2T3_mutex);
	fflush(stdout);
	//print the id of the current task
  	printf(" ]3 "); fflush(stdout);

  	blocking3 = 1000000000*(start_4.tv_sec - end_4.tv_sec) + (start_4.tv_nsec-end_4.tv_nsec);
}

void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i = 0;
  	for (i = 0; i < 100; i++)
    	{
      		task3_code();

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec = next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec = next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

void task4_code()
{
	//print the id of the current task
  	printf(" 4[ "); fflush(stdout);
  	pthread_mutex_lock(&T1T4_mutex);
  	clock_gettime(CLOCK_REALTIME, &start_5);
  	wastetime();
	printf("%d", T1T4); 
	clock_gettime(CLOCK_REALTIME, &end_5);
	pthread_mutex_unlock(&T1T4_mutex);
	fflush(stdout);
	//print the id of the current task
  	printf(" ]4 "); fflush(stdout);
  	blocking4 = 1000000000*(start_5.tv_sec - end_5.tv_sec) + (start_5.tv_nsec-end_5.tv_nsec);
}

void *task4( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i = 0;
  	for (i = 0; i < 100; i++)
    	{
      		task4_code();

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec = next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec = next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}


