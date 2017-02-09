/*------------------------------------------------------------

  Copyright (C) 1999-2016 by Sebastian Will.
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

/* ***********************************************************
 *
 * options -- an interface for getopt_long
 *
 * ***********************************************************/

#include "options.hh"
#include "aux.hh"
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include <sstream>

#include <iostream>
#include <vector>

#include "package_constants.hh"

namespace LocARNA {

    //! string holding for error message
    std::string O_error_msg;

    /* prototypes */

    /**
     * @brief decode a string according to arg_type, for internal use
     *
     * @param argument
     * @param arg_type
     * @param optarg
     *
     * @return success
     */
    bool
    decode_argument(void *argument, int arg_type, const std::string &optarg);

    /**
     * Counts the options in (null-terminated) array
     *
     * @param options array of options
     * @see process_options()
     *
     * @return number of options in array
     */
    int
    count_opts(option_def *options);

    /**
     * Write option name to string
     *
     * @param options options array
     * @param i index of option to be printed
     *
     * @return option string
     */
    std::string
    sprint_option_name(option_def *options, int i);

    /**
     * Write an option as xml in Galaxy Wrapper format to string
     *
     * @param options options array
     * @param i index of option to be printed
     *
     * @return xml option string
     */
    std::string
    sprint_option_xml(option_def *options, int i);

    /**
     * Print option name to string buffer, marking whether it is optional by []
     *
     * @param options options array
     * @param i index of option to be printed
     *
     * @return option string
     */
    std::string
    sprint_option_name_opt(option_def *options, int i);

    /**
     * \brief process options
     *
     * Processes an argument string as given to main() according to an
     * array of options descriptions of type option_def.
     *
     * @param argc argument counter from main()
     * @param argv argument vector from main()
     * @param options array describing the options
     *
     * \return TRUE   if all required options are read
     *         FALSE  if a required option is missing,
     *                or unknown options occured
     *
     * @note last entry in definitions array should contain only 0's (in
     *   fact longname==0 and shortname==0 is tested)
     *
     */
    bool
    process_options(int argc, char *argv[], option_def *options) {
        int i, k;     /* counter */
        int num_opts; /* number of options in options[] */

        int long_index; /* index in long_opts */

        char c;

        int success; // flag indicating success of decoding an option

        std::string short_opts;
        std::vector<struct option> long_opts;
        std::vector<bool>
            is_set; // boolean var for each option: is the option set?

        int index;

        num_opts = count_opts(options);

        // std::cout <<"NUM Opts: "<<num_opts<<std::endl;

        long_opts.resize(num_opts + 1);
        is_set.resize(num_opts);

        /* generate short options string and long options struct */
        for (i = 0, k = 0; i < num_opts; ++i) {
            /* short options */
            if (options[i].shortname) {
                short_opts += options[i].shortname;
                if (options[i].argument != 0) { /* with argument */
                    short_opts += ':';
                }
            }

            /* long options */
            if (options[i].longname != "") {
                long_opts[k].name =
                    const_cast<char *>(options[i].longname.c_str());
                if (options[i].argument == 0)
                    long_opts[k].has_arg = no_argument;
                else {
                    long_opts[k].has_arg = required_argument;
                }

                // in combination the following assignments control how
                // the index of a long option is returned
                long_opts[k].flag = &index; // causes getopt_long to load val
                                            // into variable index
                long_opts[k].val =
                    i; // sets val to index of long option in options array

                ++k;
            }
        }

        /* clear option flags */
        for (i = 0; i < num_opts; ++i) {
            if (options[i].flag)
                *(options[i].flag) = FALSE;
            is_set[i] = FALSE;
        }

        /* set default values */
        for (i = 0; i < num_opts; ++i) {
            if (options[i].arg_type > O_SECTION) {
                if (options[i].deflt != O_NODEFAULT) {
                    is_set[i] = TRUE;
                    success =
                        decode_argument(options[i].argument,
                                        options[i].arg_type, options[i].deflt);
                    if (!success) {
                        fprintf(stderr,
                                "INTERNAL ERROR. Option --%s: parsing of "
                                "default argument failed\n",
                                options[i].longname.c_str());
                        throw(failure(""));
                    }
                }
            }
        }

        /* main loop to process options */
        while ((c = getopt_long(argc, argv, short_opts.c_str(), &(long_opts[0]),
                                &long_index)) != EOF) {
            switch (c) {
                case '?':
                    return FALSE;
                case ':':
                    return FALSE;
                default:
                    if (c != 0) { /* short option */
                        /* find option index */
                        for (i = 0; i < num_opts && options[i].shortname != c;
                             ++i)
                            ;
                        assert(i < num_opts);
                        index = i;
                    } /* else: 'long option';
                         the index of long option is already written to variable
                         index by getopt_long,
                         due to our initialization of long_opts
                      */

                    if (options[index].flag)
                        *(options[index].flag) = TRUE;

                    // printf("decode %s
                    // %d\n",options[index].longname.c_str(),options[index].arg_type);

                    is_set[index] = TRUE;

                    if (options[index].arg_type != O_NO_ARG) {
                        // printf("decode optarg %s\n",optarg);

                        // for short options: remove leading =
                        if (c != 0 && optarg[0] == '=')
                            optarg++;

                        success = decode_argument(options[index].argument,
                                                  options[index].arg_type,
                                                  std::string(optarg));

                        if (!success) {
                            if (c == 0) { // options[index].longname != ""
                                O_error_msg =
                                    "Cannot parse argument of option --" +
                                    options[index].longname;
                            } else {
                                O_error_msg =
                                    "Cannot parse argument of option -";
                                O_error_msg += c;
                            }
                            return FALSE;
                        }
                    }
            }
        }

        /* read no option arguments */
        for (i = 0; i < num_opts; ++i)
            if (options[i].longname == "" && options[i].shortname == 0 &&
                options[i].arg_type > O_SECTION)
                if (optind < argc) {
                    if (options[i].flag)
                        *options[i].flag = TRUE;
                    is_set[i] = TRUE;
                    success =
                        decode_argument(options[i].argument,
                                        options[i].arg_type, argv[optind]);
                    if (!success) {
                        O_error_msg =
                            "Cannot parse argument no option argument.";
                        return FALSE;
                    }
                    optind++;
                }

        if (optind != argc) {
            O_error_msg = "Too many arguments.";
            return FALSE;
        }

        /* are mandatory option arguments set */
        for (i = 0; i < num_opts; ++i)
            if (options[i].argument != 0 && is_set[i] == FALSE &&
                options[i].flag == 0) {
                std::string head = "Mandatory option and/or argument missing: ";
                if (options[i].longname != "")
                    O_error_msg = head + "--" + options[i].longname;
                else if (options[i].shortname)
                    O_error_msg = head + "-" + options[i].shortname;
                else
                    O_error_msg = head + "<" +
                        ((options[i].argname != "") ? options[i].argname.c_str()
                                                    : "param") +
                        ">";
                return FALSE;
            }

        return TRUE;
    }

    const char *
    convert_arg_type(int arg_type) {
        switch (arg_type) {
            case O_ARG_STRING:
                return "text";
            case O_ARG_INT:
                return "integer";
            case O_ARG_FLOAT:
                return "float";
            case O_ARG_DOUBLE:
                return "float";
            case O_ARG_BOOL:
                return "boolean";
            default:
                return "UNKNOWN";
        }
    }

    void
    print_options(option_def options[]) {
        bool hide_options = false;

        int i;        /* counter */
        int num_opts; /* number of options in options[] */

        num_opts = count_opts(options);

        for (i = 0; i < num_opts; ++i) {
            if (options[i].arg_type > O_SECTION) {
                if (!hide_options) {
                    printf("  %-32s ", sprint_option_name(options, i).c_str());

                    if (options[i].flag != 0 && options[i].argument == 0) {
                        fputs((bool)(*options[i].flag) ? "ON" : "OFF", stdout);
                    } else {
                        if (options[i].flag == 0 || *options[i].flag) {
                            fputs("= ", stdout);
                            if (options[i].argument)
                                switch (options[i].arg_type) {
                                    case O_ARG_STRING:
                                        printf("\"%s\"",
                                               ((std::string *)options[i]
                                                    .argument)
                                                   ->c_str());
                                        break;
                                    case O_ARG_INT:
                                        printf("%d",
                                               *((int *)(options[i].argument)));
                                        break;
                                    case O_ARG_FLOAT:
                                        printf("%f",
                                               *((float *)(options[i]
                                                               .argument)));
                                        break;
                                    case O_ARG_DOUBLE:
                                        printf("%f",
                                               *((double *)(options[i]
                                                                .argument)));
                                        break;
                                    case O_ARG_BOOL:
                                        if (*((bool *)(options[i].argument)))
                                            fputs("true", stdout);
                                        else
                                            fputs("false", stdout);
                                        break;
                                    default:
                                        fputs("has unknown type", stdout);
                                }
                            else {
                                fputs("ON", stdout);
                            }
                        } else {
                            fputs("-", stdout);
                        }
                    }
                    fputs("\n", stdout);
                }    // end if (!hide_options)
            } else { // NEW SECTION
                hide_options = (options[i].arg_type == O_SECTION_HIDE);
                if (!hide_options) {
                    printf("%s:\n", options[i].description.c_str());
                }
            }
        }
    }

    /**
     * @brief test whether option is mandatory
     */
    bool
    mandatory(option_def *options, int i) {
        return options[i].flag == 0 && (options[i].deflt == O_NODEFAULT);
    }

    /**
     * @brief test whether option is positional
     */
    bool
    positional(option_def *options, int i) {
        return options[i].arg_type > O_SECTION &&
            options[i].arg_type != O_TEXT && !options[i].shortname &&
            options[i].longname == "";
    }

    /**
     * prints a standard usage string suited for short help output
     *
     * @param progname Name of program
     * @param options  Options array
     */
    void
    print_usage(char *progname, option_def options[], bool terse) {
        bool hide_options = false; /* true for hidden sections */
        int i;                     /* counter */
        int num_opts;              /* number of options in options[] */

        num_opts = count_opts(options);

        printf("USAGE: %s ", progname);

        if (terse) {
            fputs("[options]", stdout);
        }

        for (i = 0; i < num_opts; ++i) {
            /* options and no options*/
            if (options[i].arg_type > O_SECTION) {
                if (!hide_options && mandatory(options, i)) {
                    fputs(sprint_option_name_opt(options, i).c_str(), stdout);
                }
            } else {
                hide_options = (options[i].arg_type == O_SECTION_HIDE);
            }
        }
        fputs("\n", stdout);
    }

    void
    print_galaxy_xml(char *progname, option_def options[]) {
        int i;        /* counter */
        int num_opts; /* number of options in options[] */
        bool is_if_open = false;
        bool ignore_category = false;
        num_opts = count_opts(options);
        std::string category;

        printf("<!-- Galaxy wrapper for *%s* of the package *%s*\n-->",
               progname, PACKAGE_STRING);
        printf("<!-- Automatically generated by %s (option galaxy-xml)\n-->",
               progname);
        fputs("<!-- Please do NOT edit the generated wrapper.\n-->", stdout);
        printf("<!-- Source repository: %s\n-->\n\n", PACKAGE_VCS);

        printf(
            "<tool id=\"%s\" name=\"%s\" version=\"%s\">\n"
            "    <requirements>\n"
            "        <requirement type=\"package\" "
            "version=\"%s\">%s</requirement>\n"
            "    </requirements>\n"
            "    <stdio>\n"
            "        <exit_code range=\"1:\" />\n"
            "    </stdio>\n",
            progname, progname, VERSION, VERSION, PACKAGE_NAME);

        //======================================================================
        // print commands
        fputs("    <command><![CDATA[\n", stdout);
        printf("%s\n", progname);
        fputs(
            "        '$input1'\n"
            "        '$input2'\n"
            "        --clustal '$clustal_output'\n",
            stdout);

        for (i = 0; i < num_opts; ++i) {
            /* options and no options*/
            if (options[i].arg_type == O_SECTION) {
                if (is_if_open) {
                    fputs("        #end if\n", stdout);
                    is_if_open = false;
                }
                if (options[i].description == "cmd_only" ||
                    options[i].description.substr(0, 11) == "Input_files") {
                    ignore_category = true;
                } else // Valid parameter category
                {
                    ignore_category = false;
                    // Category name is the first token of description with
                    // space delim
                    std::string description = options[i].description;
                    size_type split_pos = description.find(' ');
                    if (split_pos == std::string::npos)
                        split_pos = description.length();
                    category = description.substr(0, split_pos);
                    assert(category.length() > 0);

                    printf("        #if $%s_selector\n", category.c_str());
                    is_if_open = true;
                }

            } else if (options[i].arg_type > O_SECTION && !ignore_category &&
                       options[i].longname != "clustal") {
                assert(category.length() > 0);
                std::string longname(options[i].longname);
                // replace all "-" in name with "_"
                while ((longname.find("-")) != std::string::npos)
                    longname.replace(longname.find("-"), 1, "_");

                printf("            --%s    $%s.%s",
                       options[i].longname.c_str(), category.c_str(),
                       longname.c_str());
                fputs("\n", stdout);
            }
        }
        if (is_if_open) {
            fputs("        #end if\n", stdout);
        }
        fputs("]]></command>\n", stdout);

        //======================================================================
        // print parameters
        ignore_category = false;
        bool is_conditional_open = false;

        fputs("<inputs>\n", stdout);

        // Find input files and write respective fields at the top
        for (i = 0; i < num_opts; ++i) {
            if (options[i].description.substr(0, 10) == "Input file") {
                printf("    %s ", sprint_option_xml(options, i).c_str());
                fputs("\n", stdout);
            }
        }

        for (i = 0; i < num_opts; ++i) {
            /* options and no options*/
            if (options[i].arg_type == O_SECTION ||
                options[i].arg_type == O_SECTION_HIDE) {
                if (is_conditional_open) {
                    fputs("    </section>\n", stdout);
                    is_conditional_open = false;
                }
                // Hidden section, Input files and cmd_only not needed condition
                if (options[i].arg_type == O_SECTION_HIDE ||
                    options[i].description == "cmd_only" ||
                    options[i].description.substr(0, 11) == "Input_files") {
                    ignore_category = true;
                } else if (options[i].description !=
                           "") // Valid parameter category
                {
                    ignore_category = false;
                    // Category name is the first token of description with
                    // space delim
                    std::string description = options[i].description;
                    size_type split_pos = description.find(' ');
                    if (split_pos == std::string::npos)
                        split_pos = description.length();
                    std::string category = description.substr(0, split_pos);
                    assert(category.length() > 0);

                    printf(
                        "    <section name=\"%s\" title=\"%s\" "
                        "expanded=\"False\">\n",
                        category.c_str(), category.c_str());
                    is_conditional_open = true;
                }

            } else if (options[i].arg_type > O_SECTION && !ignore_category &&
                       options[i].longname != "clustal") {
                printf("        %s ", sprint_option_xml(options, i).c_str());
                fputs("\n", stdout);
            }
        }

        if (is_conditional_open) {
            fputs("    </section>\n", stdout);
        }
        fputs(
            "    </inputs>\n"
            "    <outputs>\n"
            "        <data format=\"clustal\" name=\"clustal_output\" "
            "label=\"CLUSTAL outfile \"/>\n"
            "    </outputs>\n",
            stdout);

        //======================================================================
        // print tests

        fputs(
            "    <tests>\n"
            "    <!-- *******PUT TESTS HERE******** -->\n"
            "    </tests>\n"
            "\n",
            stdout);

        //======================================================================
        // print help
        fputs("    <help><![CDATA[\n", stdout);
        printf("This tool is part of the package %s - %s.\n\n", PACKAGE_NAME,
               PACKAGE_SHORT_DESCRIPTION);
        printf("%s\n\n", PACKAGE_LONG_DESCRIPTION);
        printf(
            "Please find the online documentation at\n"
            ".. __: %s\n\n",
            PACKAGE_URL);
        printf(
            "The software package is available for download at\n"
            ".. __: %s\n",
            PACKAGE_URL);
        fputs("    ]]></help>\n\n", stdout);

        fputs(
            "    <citations>\n"
            "        <citation "
            "type=\"doi\">10.1371/journal.pcbi.0030065</citation>\n"
            "        <citation "
            "type=\"doi\">10.1093/bioinformatics/btv185</citation>\n"
            "        <citation "
            "type=\"doi\">10.1186/s12859-014-0404-0</citation>\n"
            "    </citations>\n"
            "</tool>\n",
            stdout);
    }

    void
    print_wrapped(std::string s, size_t offset, size_t width) {
        size_t tolerance = 16;

        std::string t;
        if (offset + s.length() > width) {
            t = s.substr(0, width - offset);
            s = s.substr(width - offset);

            bool goodsep = false;

            if (s[0] != ' ') {
                size_t pos = t.rfind(' ');
                if (pos + 1 == t.length()) {
                    goodsep = true;
                } else if (pos + tolerance > t.length()) {
                    s = t.substr(pos + 1) + s;
                    t = t.substr(0, pos);
                    goodsep = true;
                }
            } else {
                s = s.substr(1);
                goodsep = true;
            }

            if (!goodsep) {
                t = t + "\\";
                s = "" + s;
            }
        } else {
            t = s;
            s = "";
        }
        fputs(t.c_str(), stdout);

        if (s.length() > 0) {
            fputs("\n", stdout);
            fputs(std::string(offset, ' ').c_str(), stdout);
            print_wrapped(s, offset, width);
        }
        return;
    }

    /**
     * prints a standard help string suited for short help output
     *
     * @param progname Name of program
     * @param options  Options array
     */
    void
    print_help(char *progname, option_def options[]) {
        bool hide_options = false;

        int i;
        int num_opts = count_opts(options);

        print_usage(progname, options, true);

        fputs("\nOptions:\n", stdout);

        for (i = 0; i < num_opts; i++) {
            if (options[i].arg_type > O_SECTION) {
                if (!hide_options && !positional(options, i)) {
                    size_t offset = 4;
                    size_t width = 77;
                    if (options[i].arg_type != O_TEXT) {
                        printf("  %-5s",
                               sprint_option_name(options, i).c_str());
                        fputs("\n", stdout);
                        fputs(std::string(offset, ' ').c_str(), stdout);
                    } else {
                        // fputs("\n",stdout);
                        offset = 2;
                        fputs(std::string(offset, ' ').c_str(), stdout);
                    }
                    if (options[i].description != "") {
                        print_wrapped(options[i].description, offset, width);
                    }
                    fputs("\n\n", stdout);
                }
            } else { // NEW SECTION
                hide_options = (options[i].arg_type == O_SECTION_HIDE);

                if (!hide_options) {
                    if (options[i].description != "cmd_only") {
                        printf("\n%s:\n", options[i].description.c_str());
                    }
                }
            }
        }
        fputs("\n", stdout);
    }

    /************************************************************/

    std::string
    sprint_option_name(option_def *options, int i) {
        std::ostringstream s;

        if (options[i].shortname)
            s << "-" << options[i].shortname;
        if (options[i].shortname && (options[i].longname != ""))
            s << ",";
        if (options[i].longname != "")
            s << "--" << options[i].longname;

        if (options[i].argument) {
            if (options[i].longname != "")
                s << "=";
            s << "<" << ((options[i].argname != "") ? options[i].argname.c_str()
                                                    : "param")
              << ">";

            if (options[i].deflt != O_NODEFAULT) {
                s << "(" << options[i].deflt.c_str() << ")";
            }
        }
        return s.str();
    }

    std::string
    sprint_option_xml(option_def *options, int i) {
        std::ostringstream s;

        s << "<param ";

        if (options[i].longname == "" &&
            options[i].argname.substr(0, 5) == "input") // Input files
        {
            s << "name=\"" << options[i].argname << "\" ";
            s << "type=\"data\" format=\"fasta,clustal\" "; // TODO: support pp
                                                            // 2.0
        } else { // Other params
            if (options[i].longname != "") {
                std::string longname(options[i].longname);

                // replace all "-" in name with "_"
                while ((longname.find("-")) != std::string::npos)
                    longname.replace(longname.find("-"), 1, "_");
                s << "name=\"" << longname << "\" ";
            }
            if (options[i].argument) {
                s << "type=\"" << convert_arg_type(options[i].arg_type)
                  << "\" ";

                if (options[i].deflt == O_NODEFAULT) {
                    s << "optional=\"True\" ";
                } else {
                    s << "value=\"" << options[i].deflt << "\" ";
                }
            } else { // O_NO_ARG
                s << "type=\"boolean\" checked=\"false\" truevalue=\"--"
                  << options[i].longname << "\" falsevalue=\"\" ";
            }
        }

        if (options[i].description != "")
            s << "label=\"" << options[i].description << "\" ";

        s << "/>";

        return s.str();
    }

    std::string
    sprint_option_name_opt(option_def *options, int i) {
        std::ostringstream s;

        s << " ";

        if (!mandatory(options, i))
            s << "[";

        if (options[i].shortname)
            s << "-" << options[i].shortname;
        if (options[i].shortname && (options[i].longname != ""))
            s << ",";
        if (options[i].longname != "")
            s << "--" << options[i].longname.c_str();

        if (options[i].argument) {
            if (options[i].longname != "")
                s << "=";
            s << "<" << ((options[i].argname != "") ? options[i].argname.c_str()
                                                    : "param")
              << ">";
        }

        if (!mandatory(options, i))
            s << "]";

        return s.str();
    }

    /**
     * \brief decode an argument according to specified type
     *
     * decode the string optarg according to the type of
     * argument arg_type and write the result to argument
     *
     * @return whether argument could be decoded
     */
    bool
    decode_argument(void *argument, int arg_type, const std::string &optarg) {
        if (argument == 0) {
            fprintf(stderr, "process_options: no argument variable\n");
            throw failure("");
        }

        int success = 0;

        switch (arg_type) {
            case O_ARG_STRING:
                *((std::string *)argument) = optarg;
                success = 1;
                break;
            case O_ARG_INT:
                // use large width limits for scanf (e.g., %20d limits ints to
                // at most 20 characters)
                success = sscanf(optarg.c_str(), "%20d", (int *)argument);
                break;
            case O_ARG_FLOAT:
                success = sscanf(optarg.c_str(), "%40f", (float *)argument);
                break;
            case O_ARG_DOUBLE:
                success = sscanf(optarg.c_str(), "%80lf", (double *)argument);
                break;
            case O_ARG_BOOL:
                *((bool *)argument) = false;

                if ((optarg == "f") || (optarg == "0") || (optarg == "false") ||
                    (optarg == "off")) {
                    success = 1;
                }

                if ((optarg == "t") || (optarg == "1") || (optarg == "true") ||
                    (optarg == "on")) {
                    *((bool *)argument) = true;
                    success = 1;
                }
                break;
            default:
                fprintf(stderr, "process_options: unknown argument type\n");
                throw failure("");
                break;
        }

        return (success == 1);
    }

    int
    count_opts(option_def *options) {
        int i;
        for (i = 0;
             !(options[i].arg_type != O_TEXT && options[i].argument == NULL &&
               options[i].flag == NULL && options[i].arg_type >= 0);
             ++i)
            ;
        return i;
    }

} // end namespace LocARNA
