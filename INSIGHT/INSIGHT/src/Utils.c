/** 
   \file Utils.c 
   Collection of general utility functions
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>



/** freeCheck
    frees an array if not null.
    @param array pointer to space to be freed
    @return 0 if array was not null, -1 otherwise
*/
int freeCheck(void* array) {
	if(array != NULL) {
		free(array);
		return 0;
	}
	
	return -1;
}
/** end of freeCheck **/



static time_t startTimeStamp;

/** startTime
    initializes global start time variable
*/
void startTime (void) {
  startTimeStamp=time(NULL);
}
/** end of startTime **/



/** printTime
    prints elapsed time since start.
    @param timestr empty string to which time is written
    @return pointer to timestr
*/
char* printTime (char* timestr) {
  time_t t;
  int h, m, s;

  t=time(NULL)-startTimeStamp;
  h=t/3600; m=(t%3600)/60; s=t-(t/60)*60;
  if(h)  sprintf(timestr,"%dh%02dm%02ds", h,m,s);
  else   sprintf(timestr,"%2dm%02ds", m,s);
  return(timestr);
}
/** end of printTime **/



/******************************************************************************************************/
/******                                          END OF FILE                                     ******/
/******************************************************************************************************/
