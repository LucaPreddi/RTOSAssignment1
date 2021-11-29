//
//	Assignment 1 of Real-Time Operating Systems class of ROBOTICS ENGINEERING first year.
//	Luca Predieri mat. 4667708
//
//	I decided to get into the problem by approaching the request first on paper, reasoni-	
//	ng with all the theory discussed in the lectures and then putting what I got on the
//	code. First of all I decided to decide where were the critical zones. We can say that
//	we have 6 different zones:
//	- Z11: writing T1T2 in memory.						
//	- Z12: writing T1T4 in memory.						
//	- Z22: writing T2T3 in memory.						
//	- Z21: reading T1T2 from the memory.						
//	- Z31: reading T2T3 from the memory.							
//	- Z41: reading T1T4 from the memory.					
//	After knowing the critical zones we can finally build β*i groups, these are:	
//	- β*1 = {Z21, Z41}							
//	- β*2 = {Z31, Z41}							
//	- β*3 = {Z41}								
//	- β*4 = {0[empty]}								
//	Then we can get from these groups the blocking time for each single task by taking the 
//	highest value inside the group. We need the blocking time to calculate the U for each 
//	single task and by computing even the Ulub we can get if the tasks are schedulable. we
//	need to do this because we're using priority ceiling with the semaphores, so we need
//	to use this kind of formula. After calculating the schedulability we can assign the pr-
//	iority to semaphores, as long as we have 3 critical zones (3 different memory location)
//	we have three different mutexes:						
//	- PT1T2: with the priority of the task 1.				
//	- PT1T4: with the priority of the task 1.					
//	- PT2T3: with the priority of the task 2.			
//	then we can start the tasks. What we can see is a the wastingtime() function that I ma-
//	de to lose some time, so we get a longer computational time (U↑) in order to see if th-
//	e check of the schedulability actually works.					

// Compile with: 
// g++ assignment.cpp -pthread -o assignment

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>

// Initialisations of different parameters, like number of loops in wastetime()
// and number of periodic and aperiodic tasks.

#define INNERLOOP 100
#define OUTERLOOP 100

#define NPERIODICTASKS 4
#define NAPERIODICTASKS 0
#define NTASKS NPERIODICTASKS + NAPERIODICTASKS

//Declaring functions of periodic tasks.

void task1_code();
void task2_code();
void task3_code();
void task4_code();

// Declaring U and Ulub.
 
double U, Ulub;

// Characteristic function of the thread.

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

// Declaring time structs for studying missed deadlines and critical zones of the
// tasks.

struct timespec start_1, end_1, start_2, end_2, start_3, end_3, start_4, end_4, start_5, end_5, start_6, end_6;
double Z11, Z12, Z22, Z21, Z31, Z41;
double B1 = 0, B2 = 0, B3 = 0, B4 = 0;
double betastar1[2], 
	   betastar2[2], 
	   betastar3[1], 
	   betastar4[0];

// Declaring different arrays for the code.

long int periods[NTASKS];					// Periods of the tasks.
struct timespec next_arrival_time[NTASKS];	// Next arrival times of tasks.
double WCET[NTASKS];						// Worst case execution time of tasks.
pthread_attr_t attributes[NTASKS];			// Attributes of the threads.
pthread_t thread_id[NTASKS];				// IDs of threads.
struct sched_param parameters[NTASKS];		// Parameters of the threads.
int missed_deadlines[NTASKS];				// Missed deadlines.

// Building wastetime() function which will be used to waste some time during
// the execution of the tasks, this function can easily change schedulability
// by changing the parameters of INNERLOOP and OUTERLOOP.

void wastetime(){
	float guardie;
	for (int i =0; i<INNERLOOP; i++){
		for (int j = 0; j < OUTERLOOP; j++){
				guardie = rand()*rand();		// Just an instruction to waste time.
		}
	}
}



int main()
{
  	// Set different task periods in nanoseconds.
	// The first task has period 200 millisecond.
	// The second task has period 160 millisecond.
	// The third task has period 100 millisecond.
	// The fourth task has period 80 millisecond.

  	periods[3] = 200000000;	// In nanoseconds.
  	periods[2] = 160000000;	// In nanoseconds.
  	periods[1] = 100000000;	// In nanoseconds.
  	periods[0] = 80000000;	// In nanoseconds.

	//Assigning a name to the maximum and the minimum priority in the system.

  	struct sched_param priomax;
  	priomax.sched_priority = sched_get_priority_max(SCHED_FIFO);
  	struct sched_param priomin;
  	priomin.sched_priority = sched_get_priority_min(SCHED_FIFO);

	// Checking that the main thread is executed with superuser privileges
	// before doing anything else, also setting the maximum priority to the 
	// current thread. All this process need to be done in superuser mode.

  	if (getuid() == 0) pthread_setschedparam(pthread_self(), SCHED_FIFO, &priomax);

  	// Execute all tasks each one for each single time and repeat them 50 times in order
  	// to find the real worst case execution time.

  	for (int i = 0; i<NTASKS; i++){

			struct timespec time_1, time_2;

	 	    if (i==0) {
	 	    	for (int j = 0; j<50; j++){
		 	    	clock_gettime(CLOCK_REALTIME, &time_1);
		 			task1_code();
		 	    	clock_gettime(CLOCK_REALTIME, &time_2);
		 	    	double realtime = 1000000000*(time_2.tv_sec - time_1.tv_sec)
					      +(time_2.tv_nsec-time_1.tv_nsec);
					if (realtime > WCET[i]){
						WCET[i] = realtime;
					}
				}
	 	    }
	        if (i==1) {
	        	for (int j = 0; j<50; j++){
		 	    	clock_gettime(CLOCK_REALTIME, &time_1);
		 			task2_code();
		 	    	clock_gettime(CLOCK_REALTIME, &time_2);
		 	    	double realtime = 1000000000*(time_2.tv_sec - time_1.tv_sec)
					      +(time_2.tv_nsec-time_1.tv_nsec);
					if (realtime > WCET[i]){
						WCET[i] = realtime;
					}
				}
	        }
	      	if (i==2) {
	      		for (int j = 0; j<50; j++){
		 	    	clock_gettime(CLOCK_REALTIME, &time_1);
		 			task3_code();
		 	    	clock_gettime(CLOCK_REALTIME, &time_2);
		 	    	double realtime = 1000000000*(time_2.tv_sec - time_1.tv_sec)
					      +(time_2.tv_nsec-time_1.tv_nsec);
					if (realtime > WCET[i]){
						WCET[i] = realtime;
					}
				}
			}
			if (i==3) {
				for (int j = 0; j<50; j++){
		 	    	clock_gettime(CLOCK_REALTIME, &time_1);
		 			task4_code();
		 	    	clock_gettime(CLOCK_REALTIME, &time_2);
		 	    	double realtime = 1000000000*(time_2.tv_sec - time_1.tv_sec)
					      +(time_2.tv_nsec-time_1.tv_nsec);
					if (realtime > WCET[i]){
						WCET[i] = realtime;
					}
				}
			}

    printf("\nWorst Case Execution Time %d=%f.\n", i, WCET[i]);
    }

    // Filling Beta*.

    betastar1[0] = Z21;
    betastar1[1] = Z41;
    betastar2[0] = Z31;
    betastar2[1] = Z41;
    betastar3[0] = Z41;

    // Computing blocking times.

    for(int i = 0; i <sizeof(betastar1)/sizeof(betastar1[0]); i++){
    	if(betastar1[i]>B1){
    		B1 = betastar1[i];
    	}
    }
    for(int i = 0; i <sizeof(betastar2)/sizeof(betastar2[0]); i++){
    	if(betastar1[i]>B2){
    		B2 = betastar2[i];
    	}
    }
    for(int i = 0; i <sizeof(betastar3)/sizeof(betastar3[0]); i++){
    	if(betastar3[i]>B3){
    		B3 = betastar3[i];
    	}
    }
    for(int i = 0; i <sizeof(betastar4)/sizeof(betastar4[0]); i++){
    	if(betastar4[i]>B4){
    		B4 = betastar4[i];
    	}
    }

    // Printing blocking times to check some errors.

    printf(" B1 %f, B2 %f, B3 %f, B4 %f. \n", B1, B2, B3, B4);

	// Calculating the schedulability of the tasks with priority ceiling. I decided
	// to calculate the U and Ulub by using a for loop.

	for(int i = 1 ; i<=NPERIODICTASKS; i++){
	
		double Ulub = i*(pow(2.0,(1.0/i)) -1);

		switch(i){

			case(1):
				U = WCET[i-1]/periods[i-1] + B1/(periods[i-1]);
				printf("U1: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(2):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + B2/(periods[i-1]);
				printf("U2: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(3):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + WCET[i-3]/periods[i-3] + B3/(periods[i-1]);
				printf("U3: %f, Ulub: %f\n", U, Ulub); fflush(stdout);
			break;

			case(4):
				U = WCET[i-1]/periods[i-1] + WCET[i-2]/periods[i-2] + WCET[i-3]/periods[i-3] + WCET[i-4]/periods[i-4] + B4/(periods[i-4]);
				printf("U4: %f, Ulub: %f\n", U, Ulub); fflush(stdout);

			break;
		}

		// If we won't have a schedulable task, the process will immediatly finish.

		if (U > Ulub){
      		printf("\n U=%lf Ulub=%lf Non schedulable Task Set, %d\n", U, Ulub, i);
      		return(-1);
    	}
	}
	printf("Schedulable test!\n");
	fflush(stdout);
  	sleep(5); 							// Just to let people read infos.

  	// Setting the minimum priority to the current thread. We need it because 
	// we will assign higher priorities to periodic threads we'll create.

  	if (getuid() == 0) pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  	// Setting the attributes of each task, including scheduling policy and prority.

  	for (int i = 0; i < NPERIODICTASKS; i++){

		// Initializing attribute structure of i-task.

      	pthread_attr_init(&(attributes[i]));

		// Setting attributes to tell the kernel that priorities and policies are explicitly chosen.

      	pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// Setting attributes to set the SCHED_FIFO policy.

		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

		// Setting the parameters to assign the priority (inversely proportional).

      	parameters[i].sched_priority = priomin.sched_priority + NTASKS - i;

		// Setting attributes and parameters of the current thread.

      	pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    }

    // Declaring and initializing attributes of the mutexes, setting the protocol of 
    // mymutexattr.

    pthread_mutexattr_t mymutexattr; pthread_mutexattr_init(&mymutexattr); 
	pthread_mutexattr_setprotocol(&mymutexattr, PTHREAD_PRIO_PROTECT);

	// Assigning priority to each single mutex.

	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority); 
	pthread_mutex_init(&T1T2_mutex, &mymutexattr); 
	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[1].sched_priority); 
	pthread_mutex_init(&T2T3_mutex, &mymutexattr);
	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority);
	pthread_mutex_init(&T1T4_mutex, &mymutexattr);

	// Declaring the variable containing return values of pthread_create.

  	int iret[NTASKS];

	// Declaring variables to read the current time.

	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1);

  	// Setting the next arrival time for each single task.

  	for (int i = 0; i < NPERIODICTASKS; i++){
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i];
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}

	// Creating all threads.

  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
  	iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// Joining all threads.

  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
  	pthread_join( thread_id[3], NULL);

  	// Printing all the variables.

  	printf("\nThe variables are:\n");
  	printf("T1T2 = %d\n", T1T2);
  	printf("T1T4 = %d\n", T1T4);
  	printf("T2T3 = %d\n\n", T2T3);

  	// Printing all missed deadlines.

  	for (int i = 0; i < NTASKS; i++){
      	printf ("Missed Deadlines Task %d=%d\n", i, missed_deadlines[i]);
		fflush(stdout);
    }

    // Closing the process.

  	exit(0);
}

// Task 1 application.

void task1_code()
{
	// Print the ID of the current task.

  	printf(" 1[ "); fflush(stdout);

  	// Entering in the critical zones.

	pthread_mutex_lock(&T1T2_mutex);		// Locking the mutex.
	clock_gettime(CLOCK_REALTIME, &start_1);	// Starting the calculating time.
	wastetime();
	T1T2 += 1;					// Writing T1T2.
	clock_gettime(CLOCK_REALTIME, &end_1);	// Finishing the calculating time.
	pthread_mutex_unlock(&T1T2_mutex);		// Unlocking the mutex.	

	pthread_mutex_lock(&T1T4_mutex);		// Locking the mutex.
	clock_gettime(CLOCK_REALTIME, &start_2);	// Starting the calculating time.
	wastetime();	
	T1T4 += 1;					// Writing T1T4.
	clock_gettime(CLOCK_REALTIME, &end_2);	// Finishing the calculating time.
	pthread_mutex_unlock(&T1T4_mutex);		// Unlocking the mutex.

	// Calculating the time for each critical zone.

	Z11 = 1000000000*(end_1.tv_sec - start_1.tv_sec) + (end_1.tv_nsec-start_1.tv_nsec);
	Z12 = 1000000000*(end_2.tv_sec - end_2.tv_sec) + (end_2.tv_nsec-start_2.tv_nsec);

  	// Print the ID of the current task.

  	printf(" ]1 "); fflush(stdout);
}

// Thread code for task 1.

void *task1( void *ptr)
{
	// Setting thread affinity.

	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

   	// Executing the task 100 times.

  	for (int i = 0; i < 100; i++)
    	{
      			// Executing application code.

			task1_code();

			// Checking for missing deadlines.

			struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[0].tv_sec*1000000000 + next_arrival_time[0].tv_nsec;

			if(currenttimed>arrivaltime) missed_deadlines[0]++;

			// Sleeping until the end of the current period (which is also the start of the
			// new one.

			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

			// The thread is ready again.

			long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
			next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

// Task 2 application.

void task2_code()
{
	// Print the ID of the current task.

  	printf(" 2[ "); fflush(stdout);

  	// Entering in the critical zones.

	pthread_mutex_lock(&T2T3_mutex);		// Locking the mutex.
	clock_gettime(CLOCK_REALTIME, &start_3);	// Starting the calculating time.
	wastetime();
	T2T3 += 1;					// Writing T2T3;
	clock_gettime(CLOCK_REALTIME, &end_3);	// Finishing the calculating time.
	pthread_mutex_unlock(&T2T3_mutex);		// Unlocking the mutex

	pthread_mutex_lock(&T1T2_mutex);		// Locking the mutex.
	clock_gettime(CLOCK_REALTIME, &start_4);	// Starting the calculating time.
	wastetime();
	int readingT1T2 = T1T2;						
	clock_gettime(CLOCK_REALTIME, &end_4);	// Finishing the calculating time.
	pthread_mutex_unlock(&T1T2_mutex);		// Unlocking the mutex.

	// Print the ID of the current task.

  	printf(" ]2 "); fflush(stdout);

  	// Calculating the time for each critical zone.

  	Z22 = 1000000000*(end_3.tv_sec - start_3.tv_sec) + (end_3.tv_nsec-start_3.tv_nsec);
  	Z21 = 1000000000*(end_4.tv_sec - start_4.tv_sec) + (end_4.tv_nsec-start_4.tv_nsec);
}

// Thread code for task 2.

void *task2( void *ptr )
{
	// Setting thread affinity.

	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Executing the task 100 times.

  	for (int i = 0; i < 100; i++)
    	{
    		// Executing application code.

      		task2_code();

      		// Checking for missing deadlines.

      		struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[1].tv_sec*1000000000 + next_arrival_time[1].tv_nsec;

			if(currenttimed>arrivaltime) missed_deadlines[1]++;

			// Sleeping until the end of the current period (which is also the start of the
			// new one

			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);

			long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
			next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

// Task 3 application.

void task3_code()
{
	// Print the ID of the current task.

  	printf(" 3[ "); fflush(stdout);

  	// Entering in the critical zones.

  	pthread_mutex_lock(&T2T3_mutex);		// Locking the mutex.
  	clock_gettime(CLOCK_REALTIME, &start_5);	// Starting the calculating time.
  	wastetime();
	int readingT2T3 = T2T3; 			// Reading T2T3.
	clock_gettime(CLOCK_REALTIME, &end_5);	// Finishing the calculating time.
	pthread_mutex_unlock(&T2T3_mutex);		// Unlocking the mutex.

	// Print the ID of the current task.

  	printf(" ]3 "); fflush(stdout);

	// Calculating the time for each critical zone.

  	Z31 = 1000000000*(end_5.tv_sec - start_5.tv_sec) + (end_5.tv_nsec-start_5.tv_nsec);
}

// Thread code for task 3.

void *task3( void *ptr)
{
	// Setting thread affinity.

	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Executing the task 100 times.

  	for (int i = 0; i < 100; i++)
    	{	
    		// Executing application code.

      		task3_code();

      		// Checking for missing deadlines.

      		struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[2].tv_sec*1000000000 + next_arrival_time[2].tv_nsec;

			if(currenttimed>arrivaltime) missed_deadlines[2]++;

			// Sleeping until the end of the current period (which is also the start of the
			// new one.

			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);

			// The thread is ready again.

			long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
			next_arrival_time[2].tv_nsec = next_arrival_nanoseconds%1000000000;
			next_arrival_time[2].tv_sec = next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

// Task 2 application.

void task4_code()
{

	// Print the ID of the current task.

  	printf(" 4[ "); fflush(stdout);

  	// Entering in the critical zones.

  	pthread_mutex_lock(&T1T4_mutex);		// Locking the mutex.
  	clock_gettime(CLOCK_REALTIME, &start_6);	// Starting the calculating time.
  	wastetime();
	int readingT1T4 = T1T4; 			// Reading T1T4.
	clock_gettime(CLOCK_REALTIME, &end_6);	// Finishing the calculating time.
	pthread_mutex_unlock(&T1T4_mutex);		// Unlocking the mutex.

  	// Calculating the time for each critical zone.

  	Z41 = 1000000000*(end_6.tv_sec - start_6.tv_sec) + (end_6.tv_nsec-start_6.tv_nsec);

  	// Print the ID of the current task.

  	printf(" ]4 "); fflush(stdout);
}

// Thread code for task 4.

void *task4( void *ptr)
{
	// Setting thread affinity.

	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Executing the task 100 times.

  	for (int i = 0; i < 100; i++){

  	// Executing application code.

      	task4_code();

      	// Checking for missing deadlines.

      	struct timespec currenttime;
		clock_gettime(CLOCK_REALTIME, &currenttime);
		double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
		double arrivaltime = next_arrival_time[3].tv_sec*1000000000 + next_arrival_time[3].tv_nsec;

		if(currenttimed>arrivaltime) missed_deadlines[3]++;

		// Sleeping until the end of the current period (which is also the start of the
		// new one.

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);

		// The thread is ready again.

		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec = next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec = next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}
