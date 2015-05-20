#ifndef UTILS_H
#define UTILS_H
/** 
   \file Utils.h
   Header file for collection of general utility functions
*/



#define FILENAME_LEN  200

/*--- macro for custom log function using log1p which is faster for values around 1 ---*/
#define customLog(x) ((1.0-x<1.0)?log1p(x-1):log(x))

/** freeCheck
    frees an array if not null.
    @param array pointer to space to be freed
    @return 0 if array was not null, -1 otherwise
*/
int freeCheck(void* array);



/** startTime
    initializes global start time variable
*/
void startTime (void);



/** printTime
    prints elapsed time since start.
    @param timestr empty string to which time is written
    @return pointer to timestr
*/
char* printTime (char* timestr);


#endif
