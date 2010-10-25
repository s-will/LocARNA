/*------------------------------------------------------------

Copyright (C) 1999 by Sebastian Will.
All Rights Reserved.
  
Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.
  
THE AUTHORS AND PUBLISHER MAKE NO REPRESENTATIONS OR
WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
THIS SOFTWARE OR ITS DERIVATIVES.

------------------------------------------------------------*/

/************************************************************
 *
 * options -- an interface for getopt_long
 *
 ************************************************************/

#include "options.hh"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


#include <iostream>

#define MYBUFSIZ 1024

char O_error_msg_buffer[MYBUFSIZ];
char *O_error_msg=O_error_msg_buffer;

/* prototype */
/* decode a string according to arg_type, for internal use */ 
bool decode_argument(void *argument, int arg_type, const std::string &optarg);
int count_opts(option_def *options);
char *sprint_option_name(char *buf,option_def *options,int i);
char *sprint_option_name_opt(char *buf,option_def *options,int i);

char buf[512]; /* buffer used for constructing strings */

/************************************************************
 *
 * process_options()
 * interface to process options 
 *
 * returns TRUE   if all required options are read
 *         FALSE  if a required option is missing,
 *                or unknown options occured 
 *
 * last entry in definitions array should contain only 0's
 *   (in fact longname==0 and shortname==0 is tested)
 *
 ************************************************************/
bool
process_options(int argc, char *argv[], option_def *options) {
    int i,j,k;        /* counter */
    int num_opts; /* number of options in options[] */

    int long_index; /* index in long_opts */
  
    char c;
  
    int success; // flag indicating success of decoding an option

    char *short_opts; 
    struct option *long_opts;
    bool *is_set; // boolean var for each option: is the option set?

    int index;

    num_opts=count_opts(options);

    //std::cout <<"NUM Opts: "<<num_opts<<std::endl;

    /* alloc sufficient memory */
    short_opts = (char *)malloc(num_opts*2 + 1);
    long_opts  = (struct option *)malloc((num_opts+1)*sizeof(struct option));
    is_set     = (bool *)malloc(num_opts*sizeof(bool)); 
			
    /* generate short options string and long options struct */
    for (i=0,j=0,k=0; i<num_opts; ++i) {
	/* short options */
	if (options[i].shortname) {
	    short_opts[j] = options[i].shortname;
	    ++j;
	    if (options[i].argument!=0) { /* with argument */
		short_opts[j]=':';
		++j;
	    }
	}
	/* long options */
	if (options[i].longname!="") {
	  long_opts[k].name = options[i].longname.c_str();
	    if (options[i].argument==0) 
		long_opts[k].has_arg = no_argument;
	    else {
		long_opts[k].has_arg = required_argument;
	    }
	    long_opts[k].flag = &index;
	    long_opts[k].val  = i;
	    ++k;
	}
    }
  
    /* clear option flags */
    for (i=0; i<num_opts; ++i) {
	if (options[i].flag) *(options[i].flag) = FALSE;
	is_set[i] = FALSE;
    }
    
    /* set default values */
    for (i=0; i<num_opts; ++i) {
	if (options[i].arg_type>O_SECTION) {
	    if (options[i].deflt!=O_NODEFAULT) {
		is_set[i] = TRUE;
		success=
		decode_argument(options[i].argument,
				options[i].arg_type,
				options[i].deflt);
		if (!success) {
		    printf("INTERNAL ERROR. Option --%s: parsing of default argument failed\n",options[i].longname.c_str());
		    exit(-1);
		}
	    }
	}
    }

    /* main loop to process options */
    while ((c=getopt_long(argc,argv,short_opts,long_opts,&long_index))!=EOF) {
	switch (c) {
	    case '?':
		return FALSE;
	    case ':':
		return FALSE;
	    default:
		if (c!=0) { /* short option */
		    /* find option index */
		    for (i=0; i<num_opts && options[i].shortname != c; ++i)
			;
		    assert(i < num_opts);
		    index = i;
		} /* else: long option */
		
		if (options[index].flag) *(options[index].flag)=TRUE;
      
		// printf("decode %s %d\n",options[index].longname.c_str(),options[index].arg_type);
		
		is_set[index] = TRUE;
		
		if (options[index].arg_type!=O_NO_ARG) {
		  
		  //printf("decode optarg %s\n",optarg);
		  
		  // for short options: remove leading =
		  if (c!=0 && optarg[0]=='=') optarg++;
		  
		  success=
		    decode_argument(options[index].argument,
				    options[index].arg_type,
				    std::string(optarg));
		  
		  if (!success) {
		    if (c==0) { // options[index].longname != ""
		      snprintf(O_error_msg,MYBUFSIZ,"Cannot parse argument of option --%s.\n",options[index].longname.c_str());
		    } else {
		      snprintf(O_error_msg,MYBUFSIZ,"Cannot parse argument of option -%c.\n",c);
		    }
		    return FALSE;
		  }

		}
	}
    }

    /* read no option arguments */
    for (i=0; i<num_opts; ++i)
	if (options[i].longname==""
	    && options[i].shortname==0
	    && options[i].arg_type>O_SECTION)
	    if (optind<argc) {
		if (options[i].flag) *options[i].flag=TRUE;
		is_set[i] = TRUE;
		success=
		decode_argument(options[i].argument,
				options[i].arg_type,
				argv[optind]);
		if (!success) {
		  snprintf(O_error_msg,MYBUFSIZ,"Cannot parse argument no option argument.");
		  return FALSE;
		}
		optind++;
	    }

    if (optind != argc) {
	snprintf(O_error_msg,MYBUFSIZ,"Too many arguments.\n");
	return FALSE;
    }
    
    
    /* are mandatory option arguments set */
    for (i=0; i<num_opts; ++i)
	if (options[i].argument != 0
	    && is_set[i]        == FALSE
	    && options[i].flag  == 0) 
	{
	  char *head=(char *)"Mandatory option and/or argument missing: ";
	    if (options[i].longname!="")
		snprintf(O_error_msg,MYBUFSIZ,"%s--%s\n",head,options[i].longname.c_str());
	    else if (options[i].shortname)
		snprintf(O_error_msg,MYBUFSIZ,"%s-%c\n",head,options[i].shortname);
	    else
		snprintf(O_error_msg,MYBUFSIZ,"%s<%s>\n",head,
			 (options[i].argname!="")?options[i].argname.c_str():"param");
	    return FALSE;
	}
  
    /* free allocated memory */
    free(short_opts);
    free(long_opts);
    free(is_set);
  
    return TRUE;
}


void
print_options(option_def options[]) {
    bool hide_options=false;

    int i;        /* counter */
    int num_opts; /* number of options in options[] */

    num_opts=count_opts(options);

    for (i=0; i < num_opts; ++i) {
	if (options[i].arg_type > O_SECTION) {
	    if (!hide_options) {
		printf("  %-32s ", sprint_option_name(buf,options,i));
		
		if (options[i].flag!=0
		    && options[i].argument==0)
		{
		    printf((bool)(*options[i].flag)?"ON":"OFF");
		} else {
		    if (options[i].flag==0 || *options[i].flag) {
			printf("= ");
			if (options[i].argument)
			    switch (options[i].arg_type) {
				case O_ARG_STRING:
				    printf ("\"%s\"",*((char **)options[i].argument));
				    break;
				case O_ARG_INT:
				    printf ("%d",*((int*)(options[i].argument)));
				    break;
				case O_ARG_FLOAT:
				    printf ("%f",*((float*)(options[i].argument)));
				    break;
				case O_ARG_DOUBLE:
				    printf ("%f",*((double*)(options[i].argument)));
				    break;
				case O_ARG_BOOL:
				    if (*((bool*)(options[i].argument))) printf("true");
				    else printf("false");
				    break;
				default:
				    printf ("has unknown type");
			    }
			else {
			    printf("ON");
			}
		    } else { 
			printf("-");
		    }
		}
		printf("\n");
	    } // end if (!hide_options)
	} else { //NEW SECTION
	    hide_options = (options[i].arg_type == O_SECTION_HIDE);
	    if (!hide_options) {
		printf("=== %s ===\n",options[i].description.c_str());
	    }
	}
    }
}


/************************************************************
 * void
 * print_usage()
 *
 * prints a standard usage string suited for short help output 
 *
 ************************************************************/
void
print_usage(char *progname, option_def options[]) {
    bool hide_options=false; /* true for hidden sections */
    int i;        /* counter */
    int num_opts; /* number of options in options[] */

    num_opts=count_opts(options);
  
    printf("%s ", progname);
   
    for (i=0; i < num_opts; ++i) {
	/* options and no options*/
	if (options[i].arg_type>O_SECTION) {
	    if (!hide_options) {
		printf("%s",sprint_option_name_opt(buf,options,i));
	    }
	} else {
	    hide_options = (options[i].arg_type == O_SECTION_HIDE);
	    if (!hide_options) printf(" ");
	}
    }
}

void
print_help(char *progname, option_def options[]) {
    bool hide_options=false;

    int i;
    int num_opts = count_opts(options);

    printf("USAGE: "); print_usage(progname, options); 
    puts("\n");

    puts("Options:");

    for (i=0; i<num_opts; i++) {
	if (options[i].arg_type>O_SECTION) {
	    if (!hide_options) {
		printf("    %-33s ",sprint_option_name(buf,options,i)); 
		
		if (options[i].description!="")
		    printf("%s",options[i].description.c_str());
		
		printf("\n");
	    }
	} else { //NEW SECTION
	    hide_options = (options[i].arg_type == O_SECTION_HIDE);
	  	  
	    if (!hide_options) {
		puts("");
		printf(" === %s ===\n",options[i].description.c_str());
	    }
	}
    }
    puts("");
}

/************************************************************/


char *sprint_option_name(char *buf,option_def *options,int i) {
    char *start=buf;
    if (options[i].shortname) buf += sprintf(buf,"-%c",options[i].shortname);
    if (options[i].shortname && (options[i].longname!="")) buf += sprintf(buf,",");
    if (options[i].longname!="") buf += sprintf(buf,"--%s",options[i].longname.c_str());
  
    if (options[i].argument) {
	if (options[i].longname!="") buf+=sprintf(buf,"=");
	buf += sprintf(buf,"<%s>", (options[i].argname!="")?options[i].argname.c_str():"param");
	if (options[i].deflt!=O_NODEFAULT) buf += sprintf(buf,"(%s)",options[i].deflt.c_str());
    }
    return start;
}

char *sprint_option_name_opt(char *buf,option_def *options,int i) {
    char *start=buf;
    bool mandatory = options[i].flag==0 && (options[i].deflt!=O_NODEFAULT) ;
  
    buf += sprintf(buf," ");

    if (!mandatory) buf += sprintf(buf,"[");
  
    if (options[i].shortname) buf += sprintf(buf,"-%c",options[i].shortname);
    if (options[i].shortname && (options[i].longname!="")) buf += sprintf(buf,",");
    if (options[i].longname!="") buf += sprintf(buf,"--%s",options[i].longname.c_str());
      
    if (options[i].argument) {
	if (options[i].longname!="") buf+=sprintf(buf,"=");
	buf += sprintf(buf,"<%s", (options[i].argname!="")?options[i].argname.c_str():"param");
	buf += sprintf(buf,">");
    }
  
    if (!mandatory) buf += sprintf(buf,"]");
    return start;
}


// decode the string <optarg> according to the type of argument <arg_type> and write the result to <argument>
// returns whether argument could be decoded
//
bool
decode_argument(void *argument, int arg_type, const std::string &optarg) {
    if (argument == 0) {
	fprintf(stderr,"process_options: no argument variable\n");
	exit(-1);
    }
  
    int success=0;

    switch(arg_type) {
    case O_ARG_STRING:
	*((std::string *)argument) = optarg;
	success=1;
	break;
    case O_ARG_INT:
	success=sscanf(optarg.c_str(), "%d", (int *)argument);
	break;
    case O_ARG_FLOAT:
	success=sscanf(optarg.c_str(), "%f", (float *)argument);
	break;
    case O_ARG_DOUBLE:
	success=sscanf(optarg.c_str(), "%lf", (double *)argument);
	break;
    case O_ARG_BOOL:
	*((bool *)argument) = false;
	
	if ((optarg=="f") || (optarg=="0") || (optarg=="false") || (optarg=="off")) {
	  success=1;
	}
	
	if ((optarg=="t")
	    || (optarg=="1")
	    || (optarg=="true")
	    || (optarg=="on")
	    ) {
	    *((bool *)argument) = true;
	    success=1;
	}
	break;
    default:
	fprintf(stderr,"process_options: unknown argument type\n");
	exit(-1);
	break;
    }
    
    return (success==1);
}

/* count entries in options */
int count_opts(option_def *options) {
    int i;
    for (i=0; !(options[i].argument == NULL && options[i].flag == NULL && options[i].arg_type>=0)
	     ; ++i)
	;
    return i;
}
