# ltmain.sh - Provide generalized library-building support services.
# NOTE: Changing this file will not affect anything until you rerun configure.
#
# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001
# Free Software Foundation, Inc.
# Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
# 
# All rights reserved.
#  
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the University of Copenhagen nor the names of its
#   contributors may be used to endorse or promote products derived from 
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

# Check that we have a working $echo.
if test "X$1" = X--no-reexec; then
  # Discard the --no-reexec flag, and continue.
  shift
elif test "X$1" = X--fallback-echo; then
  # Avoid inline document here, it may be left over
  :
elif test "X`($echo '\t') 2>/dev/null`" = 'X\t'; then
  # Yippee, $echo works!
  :
else
  # Restart under the correct shell, and then maybe $echo will work.
  exec $SHELL "$0" --no-reexec ${1+"$@"}
fi

if test "X$1" = X--fallback-echo; then
  # used as fallback echo
  shift
  cat <<EOF
$*
EOF
  exit 0
fi

# The name of this program.
progname=`$echo "$0" | ${SED} 's%^.*/%%'`
modename="$progname"

# Constants.
PROGRAM=ltmain.sh
PACKAGE=libtool
VERSION=1.4.3
TIMESTAMP=" (1.922.2.111 2002/10/23 02:54:36)"

default_mode=
help="Try \`$progname --help' for more information."
magic="%%%MAGIC variable%%%"
mkdir="mkdir"
mv="mv -f"
rm="rm -f"

# Sed substitution that helps us do robust quoting.  It backslashifies
# metacharacters that are still active within double-quoted strings.
Xsed="${SED}"' -e 1s/^X//'
sed_quote_subst='s/\([\\`\\"$\\\\]\)/\\\1/g'
# test EBCDIC or ASCII                                                         
case `echo A|od -x` in                                                         
 *[Cc]1*) # EBCDIC based system                                                
  SP2NL="tr '\100' '\n'"                                                       
  NL2SP="tr '\r\n' '\100\100'"                                                 
  ;;                                                                           
 *) # Assume ASCII based system                                                
  SP2NL="tr '\040' '\012'"                                                     
  NL2SP="tr '\015\012' '\040\040'"                                             
  ;;                                                                           
esac                                                                           

# NLS nuisances.
# Only set LANG and LC_ALL to C if already set.
# These must not be set unconditionally because not all systems understand
# e.g. LANG=C (notably SCO).
# We save the old values to restore during execute mode.
if test "${LC_ALL+set}" = set; then
  save_LC_ALL="$LC_ALL"; LC_ALL=C; export LC_ALL
fi
if test "${LANG+set}" = set; then
  save_LANG="$LANG"; LANG=C; export LANG
fi

# Make sure IFS has a sensible default
: ${IFS=" 	"}

if test "$build_libtool_libs" != yes && test "$build_old_libs" != yes; then
  echo "$modename: not configured to build any kind of library" 1>&2
  echo "Fatal configuration error.  See the $PACKAGE docs for more information." 1>&2
  exit 1
fi

# Global variables.
mode=$default_mode
nonopt=
prev=
prevopt=
run=
show="$echo"
show_help=
execute_dlfiles=
lo2o="s/\\.lo\$/.${objext}/"
o2lo="s/\\.${objext}\$/.lo/"

# Parse our command line options once, thoroughly.
while test $# -gt 0
do
  arg="$1"
  shift

  case $arg in
  -*=*) optarg=`$echo "X$arg" | $Xsed -e 's/[-_a-zA-Z0-9]*=//'` ;;
  *) optarg= ;;
  esac

  # If the previous option needs an argument, assign it.
  if test -n "$prev"; then
    case $prev in
    execute_dlfiles)
      execute_dlfiles="$execute_dlfiles $arg"
      ;;
    *)
      eval "$prev=\$arg"
      ;;
    esac

    prev=
    prevopt=
    continue
  fi

  # Have we seen a non-optional argument yet?
  case $arg in
  --help)
    show_help=yes
    ;;

  --version)
    echo "$PROGRAM (GNU $PACKAGE) $VERSION$TIMESTAMP"
    exit 0
    ;;

  --config)
    ${SED} -e '1,/^# ### BEGIN LIBTOOL CONFIG/d' -e '/^# ### END LIBTOOL CONFIG/,$d' $0
    exit 0
    ;;

  --debug)
    echo "$progname: enabling shell trace mode"
    set -x
    ;;

  --dry-run | -n)
    run=:
    ;;

  --features)
    echo "host: $host"
    if test "$build_libtool_libs" = yes; then
      echo "enable shared libraries"
    else
      echo "disable shared libraries"
    fi
    if test "$build_old_libs" = yes; then
      echo "enable static libraries"
    else
      echo "disable static libraries"
    fi
    exit 0
    ;;

  --finish) mode="finish" ;;

  --mode) prevopt="--mode" prev=mode ;;
  --mode=*) mode="$optarg" ;;

  --preserve-dup-deps) duplicate_deps="yes" ;;

  --quiet | --silent)
    show=:
    ;;

  -dlopen)
    prevopt="-dlopen"
    prev=execute_dlfiles
    ;;

  -*)
    $echo "$modename: unrecognized option \`$arg'" 1>&2
    $echo "$help" 1>&2
    exit 1
    ;;

  *)
    nonopt="$arg"
    break
    ;;
  esac
done

if test -n "$prevopt"; then
  $echo "$modename: option \`$prevopt' requires an argument" 1>&2
  $echo "$help" 1>&2
  exit 1
fi

# If this variable is set in any of the actions, the command in it
# will be execed at the end.  This prevents here-documents from being
# left over by shells.
exec_cmd=

if test -z "$show_help"; then

  # Infer the operation mode.
  if test -z "$mode"; then
    case $nonopt in
    *cc | *++ | gcc* | *-gcc* | xlc*)
      mode=link
      for arg
      do
	case $arg in
	-c)
	   mode=compile
	   break
	   ;;
	esac
      done
      ;;
    *db | *dbx | *strace | *truss)
      mode=execute
      ;;
    *install*|cp|mv)
      mode=install
      ;;
    *rm)
      mode=uninstall
      ;;
    *)
      # If we have no mode, but dlfiles were specified, then do execute mode.
      test -n "$execute_dlfiles" && mode=execute

      # Just use the default operation mode.
      if test -z "$mode"; then
	if test -n "$nonopt"; then
	  $echo "$modename: warning: cannot infer operation mode from \`$nonopt'" 1>&2
	else
	  $echo "$modename: warning: cannot infer operation mode without MODE-ARGS" 1>&2
	fi
      fi
      ;;
    esac
  fi

  # Only execute mode is allowed to have -dlopen flags.
  if test -n "$execute_dlfiles" && test "$mode" != execute; then
    $echo "$modename: unrecognized option \`-dlopen'" 1>&2
    $echo "$help" 1>&2
    exit 1
  fi

  # Change the help message to a mode-specific one.
  generic_help="$help"
  help="Try \`$modename --help --mode=$mode' for more information."

  # These modes are in order of execution frequency so that they run quickly.
  case $mode in
  # libtool compile mode
  compile)
    modename="$modename: compile"
    # Get the compilation command and the source file.
    base_compile=
    prev=
    lastarg=
    srcfile="$nonopt"
    suppress_output=

    user_target=no
    for arg
    do
      case $prev in
      "") ;;
      xcompiler)
	# Aesthetically quote the previous argument.
	prev=
	lastarg=`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`

	case $arg in
	# Double-quote args containing other shell metacharacters.
	# Many Bourne shells cannot handle close brackets correctly
	# in scan sets, so we specify it separately.
	*[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	  arg="\"$arg\""
	  ;;
	esac

	# Add the previous argument to base_compile.
	if test -z "$base_compile"; then
	  base_compile="$lastarg"
	else
	  base_compile="$base_compile $lastarg"
	fi
	continue
	;;
      esac

      # Accept any command-line options.
      case $arg in
      -o)
	if test "$user_target" != "no"; then
	  $echo "$modename: you cannot specify \`-o' more than once" 1>&2
	  exit 1
	fi
	user_target=next
	;;

      -static)
	build_old_libs=yes
	continue
	;;

      -prefer-pic)
	pic_mode=yes
	continue
	;;

      -prefer-non-pic)
	pic_mode=no
	continue
	;;

      -Xcompiler)
	prev=xcompiler
	continue
	;;

      -Wc,*)
	args=`$echo "X$arg" | $Xsed -e "s/^-Wc,//"`
	lastarg=
	save_ifs="$IFS"; IFS=','
	for arg in $args; do
	  IFS="$save_ifs"

	  # Double-quote args containing other shell metacharacters.
	  # Many Bourne shells cannot handle close brackets correctly
	  # in scan sets, so we specify it separately.
	  case $arg in
	    *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	    arg="\"$arg\""
	    ;;
	  esac
	  lastarg="$lastarg $arg"
	done
	IFS="$save_ifs"
	lastarg=`$echo "X$lastarg" | $Xsed -e "s/^ //"`

	# Add the arguments to base_compile.
	if test -z "$base_compile"; then
	  base_compile="$lastarg"
	else
	  base_compile="$base_compile $lastarg"
	fi
	continue
	;;
      esac

      case $user_target in
      next)
	# The next one is the -o target name
	user_target=yes
	continue
	;;
      yes)
	# We got the output file
	user_target=set
	libobj="$arg"
	continue
	;;
      esac

      # Accept the current argument as the source file.
      lastarg="$srcfile"
      srcfile="$arg"

      # Aesthetically quote the previous argument.

      # Backslashify any backslashes, double quotes, and dollar signs.
      # These are the only characters that are still specially
      # interpreted inside of double-quoted scrings.
      lastarg=`$echo "X$lastarg" | $Xsed -e "$sed_quote_subst"`

      # Double-quote args containing other shell metacharacters.
      # Many Bourne shells cannot handle close brackets correctly
      # in scan sets, so we specify it separately.
      case $lastarg in
      *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	lastarg="\"$lastarg\""
	;;
      esac

      # Add the previous argument to base_compile.
      if test -z "$base_compile"; then
	base_compile="$lastarg"
      else
	base_compile="$base_compile $lastarg"
      fi
    done

    case $user_target in
    set)
      ;;
    no)
      # Get the name of the library object.
      libobj=`$echo "X$srcfile" | $Xsed -e 's%^.*/%%'`
      ;;
    *)
      $echo "$modename: you must specify a target with \`-o'" 1>&2
      exit 1
      ;;
    esac

    # Recognize several different file suffixes.
    # If the user specifies -o file.o, it is replaced with file.lo
    xform='[cCFSfmso]'
    case $libobj in
    *.ada) xform=ada ;;
    *.adb) xform=adb ;;
    *.ads) xform=ads ;;
    *.asm) xform=asm ;;
    *.c++) xform=c++ ;;
    *.cc) xform=cc ;;
    *.cpp) xform=cpp ;;
    *.cxx) xform=cxx ;;
    *.f90) xform=f90 ;;
    *.for) xform=for ;;
    esac

    libobj=`$echo "X$libobj" | $Xsed -e "s/\.$xform$/.lo/"`

    case $libobj in
    *.lo) obj=`$echo "X$libobj" | $Xsed -e "$lo2o"` ;;
    *)
      $echo "$modename: cannot determine name of library object from \`$libobj'" 1>&2
      exit 1
      ;;
    esac

    if test -z "$base_compile"; then
      $echo "$modename: you must specify a compilation command" 1>&2
      $echo "$help" 1>&2
      exit 1
    fi

    # Delete any leftover library objects.
    if test "$build_old_libs" = yes; then
      removelist="$obj $libobj"
    else
      removelist="$libobj"
    fi

    $run $rm $removelist
    trap "$run $rm $removelist; exit 1" 1 2 15

    # On Cygwin there's no "real" PIC flag so we must build both object types
    case $host_os in
    cygwin* | mingw* | pw32* | os2*)
      pic_mode=default
      ;;
    esac
    if test "$pic_mode" = no && test "$deplibs_check_method" != pass_all; then
      # non-PIC code in shared libraries is not supported
      pic_mode=default
    fi

    # Calculate the filename of the output object if compiler does
    # not support -o with -c
    if test "$compiler_c_o" = no; then
      output_obj=`$echo "X$srcfile" | $Xsed -e 's%^.*/%%' -e 's%\.[^.]*$%%'`.${objext}
      lockfile="$output_obj.lock"
      removelist="$removelist $output_obj $lockfile"
      trap "$run $rm $removelist; exit 1" 1 2 15
    else
      need_locks=no
      lockfile=
    fi

    # Lock this critical section if it is needed
    # We use this script file to make the link, it avoids creating a new file
    if test "$need_locks" = yes; then
      until $run ln "$0" "$lockfile" 2>/dev/null; do
	$show "Waiting for $lockfile to be removed"
	sleep 2
      done
    elif test "$need_locks" = warn; then
      if test -f "$lockfile"; then
	echo "\
*** ERROR, $lockfile exists and contains:
`cat $lockfile 2>/dev/null`

This indicates that another process is trying to use the same
temporary object file, and libtool could not work around it because
your compiler does not support \`-c' and \`-o' together.  If you
repeat this compilation, it may succeed, by chance, but you had better
avoid parallel builds (make -j) in this platform, or get a better
compiler."

	$run $rm $removelist
	exit 1
      fi
      echo $srcfile > "$lockfile"
    fi

    if test -n "$fix_srcfile_path"; then
      eval srcfile=\"$fix_srcfile_path\"
    fi

    # Only build a PIC object if we are building libtool libraries.
    if test "$build_libtool_libs" = yes; then
      # Without this assignment, base_compile gets emptied.
      fbsd_hideous_sh_bug=$base_compile

      if test "$pic_mode" != no; then
	# All platforms use -DPIC, to notify preprocessed assembler code.
	command="$base_compile $srcfile $pic_flag -DPIC"
      else
	# Don't build PIC code
	command="$base_compile $srcfile"
      fi
      if test "$build_old_libs" = yes; then
	lo_libobj="$libobj"
	dir=`$echo "X$libobj" | $Xsed -e 's%/[^/]*$%%'`
	if test "X$dir" = "X$libobj"; then
	  dir="$objdir"
	else
	  dir="$dir/$objdir"
	fi
	libobj="$dir/"`$echo "X$libobj" | $Xsed -e 's%^.*/%%'`

	if test -d "$dir"; then
	  $show "$rm $libobj"
	  $run $rm $libobj
	else
	  $show "$mkdir $dir"
	  $run $mkdir $dir
	  status=$?
	  if test $status -ne 0 && test ! -d $dir; then
	    exit $status
	  fi
	fi
      fi
      if test "$compiler_o_lo" = yes; then
	output_obj="$libobj"
	command="$command -o $output_obj"
      elif test "$compiler_c_o" = yes; then
	output_obj="$obj"
	command="$command -o $output_obj"
      fi

      $run $rm "$output_obj"
      $show "$command"
      if $run eval "$command"; then :
      else
	test -n "$output_obj" && $run $rm $removelist
	exit 1
      fi

      if test "$need_locks" = warn &&
	 test x"`cat $lockfile 2>/dev/null`" != x"$srcfile"; then
	echo "\
*** ERROR, $lockfile contains:
`cat $lockfile 2>/dev/null`

but it should contain:
$srcfile

This indicates that another process is trying to use the same
temporary object file, and libtool could not work around it because
your compiler does not support \`-c' and \`-o' together.  If you
repeat this compilation, it may succeed, by chance, but you had better
avoid parallel builds (make -j) in this platform, or get a better
compiler."

	$run $rm $removelist
	exit 1
      fi

      # Just move the object if needed, then go on to compile the next one
      if test x"$output_obj" != x"$libobj"; then
	$show "$mv $output_obj $libobj"
	if $run $mv $output_obj $libobj; then :
	else
	  error=$?
	  $run $rm $removelist
	  exit $error
	fi
      fi

      # If we have no pic_flag, then copy the object into place and finish.
      if (test -z "$pic_flag" || test "$pic_mode" != default) &&
	 test "$build_old_libs" = yes; then
	# Rename the .lo from within objdir to obj
	if test -f $obj; then
	  $show $rm $obj
	  $run $rm $obj
	fi

	$show "$mv $libobj $obj"
	if $run $mv $libobj $obj; then :
	else
	  error=$?
	  $run $rm $removelist
	  exit $error
	fi

	xdir=`$echo "X$obj" | $Xsed -e 's%/[^/]*$%%'`
	if test "X$xdir" = "X$obj"; then
	  xdir="."
	else
	  xdir="$xdir"
	fi
	baseobj=`$echo "X$obj" | $Xsed -e "s%.*/%%"`
	libobj=`$echo "X$baseobj" | $Xsed -e "$o2lo"`
	# Now arrange that obj and lo_libobj become the same file
	$show "(cd $xdir && $LN_S $baseobj $libobj)"
	if $run eval '(cd $xdir && $LN_S $baseobj $libobj)'; then
	  # Unlock the critical section if it was locked
	  if test "$need_locks" != no; then
	    $run $rm "$lockfile"
	  fi
	  exit 0
	else
	  error=$?
	  $run $rm $removelist
	  exit $error
	fi
      fi

      # Allow error messages only from the first compilation.
      suppress_output=' >/dev/null 2>&1'
    fi

    # Only build a position-dependent object if we build old libraries.
    if test "$build_old_libs" = yes; then
      if test "$pic_mode" != yes; then
	# Don't build PIC code
	command="$base_compile $srcfile"
      else
	# All platforms use -DPIC, to notify preprocessed assembler code.
	command="$base_compile $srcfile $pic_flag -DPIC"
      fi
      if test "$compiler_c_o" = yes; then
	command="$command -o $obj"
	output_obj="$obj"
      fi

      # Suppress compiler output if we already did a PIC compilation.
      command="$command$suppress_output"
      $run $rm "$output_obj"
      $show "$command"
      if $run eval "$command"; then :
      else
	$run $rm $removelist
	exit 1
      fi

      if test "$need_locks" = warn &&
	 test x"`cat $lockfile 2>/dev/null`" != x"$srcfile"; then
	echo "\
*** ERROR, $lockfile contains:
`cat $lockfile 2>/dev/null`

but it should contain:
$srcfile

This indicates that another process is trying to use the same
temporary object file, and libtool could not work around it because
your compiler does not support \`-c' and \`-o' together.  If you
repeat this compilation, it may succeed, by chance, but you had better
avoid parallel builds (make -j) in this platform, or get a better
compiler."

	$run $rm $removelist
	exit 1
      fi

      # Just move the object if needed
      if test x"$output_obj" != x"$obj"; then
	$show "$mv $output_obj $obj"
	if $run $mv $output_obj $obj; then :
	else
	  error=$?
	  $run $rm $removelist
	  exit $error
	fi
      fi

      # Create an invalid libtool object if no PIC, so that we do not
      # accidentally link it into a program.
      if test "$build_libtool_libs" != yes; then
	$show "echo timestamp > $libobj"
	$run eval "echo timestamp > \$libobj" || exit $?
      else
	# Move the .lo from within objdir
	$show "$mv $libobj $lo_libobj"
	if $run $mv $libobj $lo_libobj; then :
	else
	  error=$?
	  $run $rm $removelist
	  exit $error
	fi
      fi
    fi

    # Unlock the critical section if it was locked
    if test "$need_locks" != no; then
      $run $rm "$lockfile"
    fi

    exit 0
    ;;

  # libtool link mode
  link | relink)
    modename="$modename: link"
    case $host in
    *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
      # It is impossible to link a dll without this setting, and
      # we shouldn't force the makefile maintainer to figure out
      # which system we are compiling for in order to pass an extra
      # flag for every libtool invokation.
      # allow_undefined=no

      # FIXME: Unfortunately, there are problems with the above when trying
      # to make a dll which has undefined symbols, in which case not
      # even a static library is built.  For now, we need to specify
      # -no-undefined on the libtool link line when we can be certain
      # that all symbols are satisfied, otherwise we get a static library.
      allow_undefined=yes
      ;;
    *)
      allow_undefined=yes
      ;;
    esac
    libtool_args="$nonopt"
    compile_command="$nonopt"
    finalize_command="$nonopt"

    compile_rpath=
    finalize_rpath=
    compile_shlibpath=
    finalize_shlibpath=
    convenience=
    old_convenience=
    deplibs=
    old_deplibs=
    compiler_flags=
    linker_flags=
    dllsearchpath=
    lib_search_path=`pwd`

    avoid_version=no
    dlfiles=
    dlprefiles=
    dlself=no
    export_dynamic=no
    export_symbols=
    export_symbols_regex=
    generated=
    libobjs=
    ltlibs=
    module=no
    no_install=no
    objs=
    prefer_static_libs=no
    preload=no
    prev=
    prevarg=
    release=
    rpath=
    xrpath=
    perm_rpath=
    temp_rpath=
    thread_safe=no
    vinfo=

    # We need to know -static, to get the right output filenames.
    for arg
    do
      case $arg in
      -all-static | -static)
	if test "X$arg" = "X-all-static"; then
	  if test "$build_libtool_libs" = yes && test -z "$link_static_flag"; then
	    $echo "$modename: warning: complete static linking is impossible in this configuration" 1>&2
	  fi
	  if test -n "$link_static_flag"; then
	    dlopen_self=$dlopen_self_static
	  fi
	else
	  if test -z "$pic_flag" && test -n "$link_static_flag"; then
	    dlopen_self=$dlopen_self_static
	  fi
	fi
	build_libtool_libs=no
	build_old_libs=yes
	prefer_static_libs=yes
	break
	;;
      esac
    done

    # See if our shared archives depend on static archives.
    test -n "$old_archive_from_new_cmds" && build_old_libs=yes

    # Go through the arguments, transforming them on the way.
    while test $# -gt 0; do
      arg="$1"
      shift
      case $arg in
      *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	qarg=\"`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`\" ### testsuite: skip nested quoting test
	;;
      *) qarg=$arg ;;
      esac
      libtool_args="$libtool_args $qarg"

      # If the previous option needs an argument, assign it.
      if test -n "$prev"; then
	case $prev in
	output)
	  compile_command="$compile_command @OUTPUT@"
	  finalize_command="$finalize_command @OUTPUT@"
	  ;;
	esac

	case $prev in
	dlfiles|dlprefiles)
	  if test "$preload" = no; then
	    # Add the symbol object into the linking commands.
	    compile_command="$compile_command @SYMFILE@"
	    finalize_command="$finalize_command @SYMFILE@"
	    preload=yes
	  fi
	  case $arg in
	  *.la | *.lo) ;;  # We handle these cases below.
	  force)
	    if test "$dlself" = no; then
	      dlself=needless
	      export_dynamic=yes
	    fi
	    prev=
	    continue
	    ;;
	  self)
	    if test "$prev" = dlprefiles; then
	      dlself=yes
	    elif test "$prev" = dlfiles && test "$dlopen_self" != yes; then
	      dlself=yes
	    else
	      dlself=needless
	      export_dynamic=yes
	    fi
	    prev=
	    continue
	    ;;
	  *)
	    if test "$prev" = dlfiles; then
	      dlfiles="$dlfiles $arg"
	    else
	      dlprefiles="$dlprefiles $arg"
	    fi
	    prev=
	    continue
	    ;;
	  esac
	  ;;
	expsyms)
	  export_symbols="$arg"
	  if test ! -f "$arg"; then
	    $echo "$modename: symbol file \`$arg' does not exist"
	    exit 1
	  fi
	  prev=
	  continue
	  ;;
	expsyms_regex)
	  export_symbols_regex="$arg"
	  prev=
	  continue
	  ;;
	release)
	  release="-$arg"
	  prev=
	  continue
	  ;;
	rpath | xrpath)
	  # We need an absolute path.
	  case $arg in
	  [\\/]* | [A-Za-z]:[\\/]*) ;;
	  *)
	    $echo "$modename: only absolute run-paths are allowed" 1>&2
	    exit 1
	    ;;
	  esac
	  if test "$prev" = rpath; then
	    case "$rpath " in
	    *" $arg "*) ;;
	    *) rpath="$rpath $arg" ;;
	    esac
	  else
	    case "$xrpath " in
	    *" $arg "*) ;;
	    *) xrpath="$xrpath $arg" ;;
	    esac
	  fi
	  prev=
	  continue
	  ;;
	xcompiler)
	  compiler_flags="$compiler_flags $qarg"
	  prev=
	  compile_command="$compile_command $qarg"
	  finalize_command="$finalize_command $qarg"
	  continue
	  ;;
	xlinker)
	  linker_flags="$linker_flags $qarg"
	  compiler_flags="$compiler_flags $wl$qarg"
	  prev=
	  compile_command="$compile_command $wl$qarg"
	  finalize_command="$finalize_command $wl$qarg"
	  continue
	  ;;
	*)
	  eval "$prev=\"\$arg\""
	  prev=
	  continue
	  ;;
	esac
      fi # test -n $prev

      prevarg="$arg"

      case $arg in
      -all-static)
	if test -n "$link_static_flag"; then
	  compile_command="$compile_command $link_static_flag"
	  finalize_command="$finalize_command $link_static_flag"
	fi
	continue
	;;

      -allow-undefined)
	# FIXME: remove this flag sometime in the future.
	$echo "$modename: \`-allow-undefined' is deprecated because it is the default" 1>&2
	continue
	;;

      -avoid-version)
	avoid_version=yes
	continue
	;;

      -dlopen)
	prev=dlfiles
	continue
	;;

      -dlpreopen)
	prev=dlprefiles
	continue
	;;

      -export-dynamic)
	export_dynamic=yes
	continue
	;;

      -export-symbols | -export-symbols-regex)
	if test -n "$export_symbols" || test -n "$export_symbols_regex"; then
	  $echo "$modename: more than one -exported-symbols argument is not allowed"
	  exit 1
	fi
	if test "X$arg" = "X-export-symbols"; then
	  prev=expsyms
	else
	  prev=expsyms_regex
	fi
	continue
	;;

      # The native IRIX linker understands -LANG:*, -LIST:* and -LNO:*
      # so, if we see these flags be careful not to treat them like -L
      -L[A-Z][A-Z]*:*)
	case $with_gcc/$host in
	no/*-*-irix* | no/*-*-nonstopux*)
	  compile_command="$compile_command $arg"
	  finalize_command="$finalize_command $arg"
	  ;;
	esac
	continue
	;;

      -L*)
	dir=`$echo "X$arg" | $Xsed -e 's/^-L//'`
	# We need an absolute path.
	case $dir in
	[\\/]* | [A-Za-z]:[\\/]*) ;;
	*)
	  absdir=`cd "$dir" && pwd`
	  if test -z "$absdir"; then
	    $echo "$modename: cannot determine absolute directory name of \`$dir'" 1>&2
	    exit 1
	  fi
	  dir="$absdir"
	  ;;
	esac
	case "$deplibs " in
	*" -L$dir "*) ;;
	*)
	  deplibs="$deplibs -L$dir"
	  lib_search_path="$lib_search_path $dir"
	  ;;
	esac
	case $host in
	*-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
	  case :$dllsearchpath: in
	  *":$dir:"*) ;;
	  *) dllsearchpath="$dllsearchpath:$dir";;
	  esac
	  ;;
	esac
	continue
	;;

      -l*)
	if test "X$arg" = "X-lc" || test "X$arg" = "X-lm"; then
	  case $host in
	  *-*-cygwin* | *-*-pw32* | *-*-beos*)
	    # These systems don't actually have a C or math library (as such)
	    continue
	    ;;
	  *-*-mingw* | *-*-os2*)
	    # These systems don't actually have a C library (as such)
	    test "X$arg" = "X-lc" && continue
	    ;;
	  *-*-openbsd* | *-*-freebsd*)
	    # Do not include libc due to us having libc/libc_r.
	    test "X$arg" = "X-lc" && continue
	    ;;
	  esac
	 elif test "X$arg" = "X-lc_r"; then
	  case $host in
	 *-*-openbsd* | *-*-freebsd*)
	    # Do not include libc_r directly, use -pthread flag.
	    continue
	    ;;
	  esac
	fi
	deplibs="$deplibs $arg"
	continue
	;;

      -module)
	module=yes
	continue
	;;

      -no-fast-install)
	fast_install=no
	continue
	;;

      -no-install)
	case $host in
	*-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
	  # The PATH hackery in wrapper scripts is required on Windows
	  # in order for the loader to find any dlls it needs.
	  $echo "$modename: warning: \`-no-install' is ignored for $host" 1>&2
	  $echo "$modename: warning: assuming \`-no-fast-install' instead" 1>&2
	  fast_install=no
	  ;;
	*) no_install=yes ;;
	esac
	continue
	;;

      -no-undefined)
	allow_undefined=no
	continue
	;;

      -o) prev=output ;;

      -release)
	prev=release
	continue
	;;

      -rpath)
	prev=rpath
	continue
	;;

      -R)
	prev=xrpath
	continue
	;;

      -R*)
	dir=`$echo "X$arg" | $Xsed -e 's/^-R//'`
	# We need an absolute path.
	case $dir in
	[\\/]* | [A-Za-z]:[\\/]*) ;;
	*)
	  $echo "$modename: only absolute run-paths are allowed" 1>&2
	  exit 1
	  ;;
	esac
	case "$xrpath " in
	*" $dir "*) ;;
	*) xrpath="$xrpath $dir" ;;
	esac
	continue
	;;

      -static)
	# The effects of -static are defined in a previous loop.
	# We used to do the same as -all-static on platforms that
	# didn't have a PIC flag, but the assumption that the effects
	# would be equivalent was wrong.  It would break on at least
	# Digital Unix and AIX.
	continue
	;;

      -thread-safe)
	thread_safe=yes
	continue
	;;

      -version-info)
	prev=vinfo
	continue
	;;

      -Wc,*)
	args=`$echo "X$arg" | $Xsed -e "$sed_quote_subst" -e 's/^-Wc,//'`
	arg=
	save_ifs="$IFS"; IFS=','
	for flag in $args; do
	  IFS="$save_ifs"
	  case $flag in
	    *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	    flag="\"$flag\""
	    ;;
	  esac
	  arg="$arg $wl$flag"
	  compiler_flags="$compiler_flags $flag"
	done
	IFS="$save_ifs"
	arg=`$echo "X$arg" | $Xsed -e "s/^ //"`
	;;

      -Wl,*)
	args=`$echo "X$arg" | $Xsed -e "$sed_quote_subst" -e 's/^-Wl,//'`
	arg=
	save_ifs="$IFS"; IFS=','
	for flag in $args; do
	  IFS="$save_ifs"
	  case $flag in
	    *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	    flag="\"$flag\""
	    ;;
	  esac
	  arg="$arg $wl$flag"
	  compiler_flags="$compiler_flags $wl$flag"
	  linker_flags="$linker_flags $flag"
	done
	IFS="$save_ifs"
	arg=`$echo "X$arg" | $Xsed -e "s/^ //"`
	;;

      -Xcompiler)
	prev=xcompiler
	continue
	;;

      -Xlinker)
	prev=xlinker
	continue
	;;

      # Some other compiler flag.
      -* | +*)
	# Unknown arguments in both finalize_command and compile_command need
	# to be aesthetically quoted because they are evaled later.
	arg=`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`
	case $arg in
	*[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	  arg="\"$arg\""
	  ;;
	esac
	;;

      *.lo | *.$objext)
	# A library or standard object.
	if test "$prev" = dlfiles; then
	  # This file was specified with -dlopen.
	  if test "$build_libtool_libs" = yes && test "$dlopen_support" = yes; then
	    dlfiles="$dlfiles $arg"
	    prev=
	    continue
	  else
	    # If libtool objects are unsupported, then we need to preload.
	    prev=dlprefiles
	  fi
	fi

	if test "$prev" = dlprefiles; then
	  # Preload the old-style object.
	  dlprefiles="$dlprefiles "`$echo "X$arg" | $Xsed -e "$lo2o"`
	  prev=
	else
	  case $arg in
	  *.lo) libobjs="$libobjs $arg" ;;
	  *) objs="$objs $arg" ;;
	  esac
	fi
	;;

      *.$libext)
	# An archive.
	deplibs="$deplibs $arg"
	old_deplibs="$old_deplibs $arg"
	continue
	;;

      *.la)
	# A libtool-controlled library.

	if test "$prev" = dlfiles; then
	  # This library was specified with -dlopen.
	  dlfiles="$dlfiles $arg"
	  prev=
	elif test "$prev" = dlprefiles; then
	  # The library was specified with -dlpreopen.
	  dlprefiles="$dlprefiles $arg"
	  prev=
	else
	  deplibs="$deplibs $arg"
	fi
	continue
	;;

      # Some other compiler argument.
      *)
	# Unknown arguments in both finalize_command and compile_command need
	# to be aesthetically quoted because they are evaled later.
	arg=`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`
	case $arg in
	*[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*|"")
	  arg="\"$arg\""
	  ;;
	esac
	;;
      esac # arg

      # Now actually substitute the argument into the commands.
      if test -n "$arg"; then
	compile_command="$compile_command $arg"
	finalize_command="$finalize_command $arg"
      fi
    done # argument parsing loop

    if test -n "$prev"; then
      $echo "$modename: the \`$prevarg' option requires an argument" 1>&2
      $echo "$help" 1>&2
      exit 1
    fi

    if test "$export_dynamic" = yes && test -n "$export_dynamic_flag_spec"; then
      eval arg=\"$export_dynamic_flag_spec\"
      compile_command="$compile_command $arg"
      finalize_command="$finalize_command $arg"
    fi

    # calculate the name of the file, without its directory
    outputname=`$echo "X$output" | $Xsed -e 's%^.*/%%'`
    libobjs_save="$libobjs"

    if test -n "$shlibpath_var"; then
      # get the directories listed in $shlibpath_var
      eval shlib_search_path=\`\$echo \"X\${$shlibpath_var}\" \| \$Xsed -e \'s/:/ /g\'\`
    else
      shlib_search_path=
    fi
    eval sys_lib_search_path=\"$sys_lib_search_path_spec\"
    eval sys_lib_dlsearch_path=\"$sys_lib_dlsearch_path_spec\"

    output_objdir=`$echo "X$output" | $Xsed -e 's%/[^/]*$%%'`
    if test "X$output_objdir" = "X$output"; then
      output_objdir="$objdir"
    else
      output_objdir="$output_objdir/$objdir"
    fi
    # Create the object directory.
    if test ! -d "$output_objdir"; then
      $show "$mkdir $output_objdir"
      $run $mkdir "$output_objdir"
      status=$?
      if test $status -ne 0 && test ! -d "$output_objdir"; then
	exit $status
      fi
    fi

    # Determine the type of output
    case $output in
    "")
      $echo "$modename: you must specify an output file" 1>&2
      $echo "$help" 1>&2
      exit 1
      ;;
    *.$libext) linkmode=oldlib ;;
    *.lo | *.$objext) linkmode=obj ;;
    *.la) linkmode=lib ;;
    *) linkmode=prog ;; # Anything else should be a program.
    esac

    specialdeplibs=
    libs=
    # Find all interdependent deplibs by searching for libraries
    # that are linked more than once (e.g. -la -lb -la)
    for deplib in $deplibs; do
      if test "X$duplicate_deps" = "Xyes" ; then
	case "$libs " in
	*" $deplib "*) specialdeplibs="$specialdeplibs $deplib" ;;
	esac
      fi
      libs="$libs $deplib"
    done
    deplibs=
    newdependency_libs=
    newlib_search_path=
    need_relink=no # whether we're linking any uninstalled libtool libraries
    notinst_deplibs= # not-installed libtool libraries
    notinst_path= # paths that contain not-installed libtool libraries
    case $linkmode in
    lib)
	passes="conv link"
	for file in $dlfiles $dlprefiles; do
	  case $file in
	  *.la) ;;
	  *)
	    $echo "$modename: libraries can \`-dlopen' only libtool libraries: $file" 1>&2
	    exit 1
	    ;;
	  esac
	done
	;;
    prog)
	compile_deplibs=
	finalize_deplibs=
	alldeplibs=no
	newdlfiles=
	newdlprefiles=
	passes="conv scan dlopen dlpreopen link"
	;;
    *)  passes="conv"
	;;
    esac
    for pass in $passes; do
      if test $linkmode = prog; then
	# Determine which files to process
	case $pass in
	dlopen)
	  libs="$dlfiles"
	  save_deplibs="$deplibs" # Collect dlpreopened libraries
	  deplibs=
	  ;;
	dlpreopen) libs="$dlprefiles" ;;
	link) libs="$deplibs %DEPLIBS% $dependency_libs" ;;
	esac
      fi
      for deplib in $libs; do
	lib=
	found=no
	case $deplib in
	-l*)
	  if test $linkmode = oldlib && test $linkmode = obj; then
	    $echo "$modename: warning: \`-l' is ignored for archives/objects: $deplib" 1>&2
	    continue
	  fi
	  if test $pass = conv; then
	    deplibs="$deplib $deplibs"
	    continue
	  fi
	  name=`$echo "X$deplib" | $Xsed -e 's/^-l//'`
	  for searchdir in $newlib_search_path $lib_search_path $sys_lib_search_path $shlib_search_path; do
	    # Search the libtool library
	    lib="$searchdir/lib${name}.la"
	    if test -f "$lib"; then
	      found=yes
	      break
	    fi
	  done
	  if test "$found" != yes; then
	    # deplib doesn't seem to be a libtool library
	    if test "$linkmode,$pass" = "prog,link"; then
	      compile_deplibs="$deplib $compile_deplibs"
	      finalize_deplibs="$deplib $finalize_deplibs"
	    else
	      deplibs="$deplib $deplibs"
	      test $linkmode = lib && newdependency_libs="$deplib $newdependency_libs"
	    fi
	    continue
	  fi
	  ;; # -l
	-L*)
	  case $linkmode in
	  lib)
	    deplibs="$deplib $deplibs"
	    test $pass = conv && continue
	    newdependency_libs="$deplib $newdependency_libs"
	    newlib_search_path="$newlib_search_path "`$echo "X$deplib" | $Xsed -e 's/^-L//'`
	    ;;
	  prog)
	    if test $pass = conv; then
	      deplibs="$deplib $deplibs"
	      continue
	    fi
	    if test $pass = scan; then
	      deplibs="$deplib $deplibs"
	      newlib_search_path="$newlib_search_path "`$echo "X$deplib" | $Xsed -e 's/^-L//'`
	    else
	      compile_deplibs="$deplib $compile_deplibs"
	      finalize_deplibs="$deplib $finalize_deplibs"
	    fi
	    ;;
	  *)
	    $echo "$modename: warning: \`-L' is ignored for archives/objects: $deplib" 1>&2
	    ;;
	  esac # linkmode
	  continue
	  ;; # -L
	-R*)
	  if test $pass = link; then
	    dir=`$echo "X$deplib" | $Xsed -e 's/^-R//'`
	    # Make sure the xrpath contains only unique directories.
	    case "$xrpath " in
	    *" $dir "*) ;;
	    *) xrpath="$xrpath $dir" ;;
	    esac
	  fi
	  deplibs="$deplib $deplibs"
	  continue
	  ;;
	*.la) lib="$deplib" ;;
	*.$libext)
	  if test $pass = conv; then
	    deplibs="$deplib $deplibs"
	    continue
	  fi
	  case $linkmode in
	  lib)
	    if test "$deplibs_check_method" != pass_all; then
	      echo
	      echo "*** Warning: Trying to link with static lib archive $deplib."
	      echo "*** I have the capability to make that library automatically link in when"
	      echo "*** you link to this library.  But I can only do this if you have a"
	      echo "*** shared version of the library, which you do not appear to have"
	      echo "*** because the file extensions .$libext of this argument makes me believe"
	      echo "*** that it is just a static archive that I should not used here."
	    else
	      echo
	      echo "*** Warning: Linking the shared library $output against the"
	      echo "*** static library $deplib is not portable!"
	      deplibs="$deplib $deplibs"
	    fi
	    continue
	    ;;
	  prog)
	    if test $pass != link; then
	      deplibs="$deplib $deplibs"
	    else
	      compile_deplibs="$deplib $compile_deplibs"
	      finalize_deplibs="$deplib $finalize_deplibs"
	    fi
	    continue
	    ;;
	  esac # linkmode
	  ;; # *.$libext
	*.lo | *.$objext)
	  if test $pass = dlpreopen || test "$dlopen_support" != yes || test "$build_libtool_libs" = no; then
	    # If there is no dlopen support or we're linking statically,
	    # we need to preload.
	    newdlprefiles="$newdlprefiles $deplib"
	    compile_deplibs="$deplib $compile_deplibs"
	    finalize_deplibs="$deplib $finalize_deplibs"
	  else
	    newdlfiles="$newdlfiles $deplib"
	  fi
	  continue
	  ;;
	%DEPLIBS%)
	  alldeplibs=yes
	  continue
	  ;;
	esac # case $deplib
	if test $found = yes || test -f "$lib"; then :
	else
	  $echo "$modename: cannot find the library \`$lib'" 1>&2
	  exit 1
	fi

	# Check to see that this really is a libtool archive.
	if (${SED} -e '2q' $lib | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then :
	else
	  $echo "$modename: \`$lib' is not a valid libtool archive" 1>&2
	  exit 1
	fi

	ladir=`$echo "X$lib" | $Xsed -e 's%/[^/]*$%%'`
	test "X$ladir" = "X$lib" && ladir="."

	dlname=
	dlopen=
	dlpreopen=
	libdir=
	library_names=
	old_library=
	# If the library was installed with an old release of libtool,
	# it will not redefine variable installed.
	installed=yes

	# Read the .la file
	case $lib in
	*/* | *\\*) . $lib ;;
	*) . ./$lib ;;
	esac

	if test "$linkmode,$pass" = "lib,link" ||
	   test "$linkmode,$pass" = "prog,scan" ||
	   { test $linkmode != lib && test $linkmode != prog; }; then
	   # Add dl[pre]opened files of deplib
	  test -n "$dlopen" && dlfiles="$dlfiles $dlopen"
	  test -n "$dlpreopen" && dlprefiles="$dlprefiles $dlpreopen"
	fi

	if test $pass = conv; then
	  # Only check for convenience libraries
	  deplibs="$lib $deplibs"
	  if test -z "$libdir"; then
	    if test -z "$old_library"; then
	      $echo "$modename: cannot find name of link library for \`$lib'" 1>&2
	      exit 1
	    fi
	    # It is a libtool convenience library, so add in its objects.
	    convenience="$convenience $ladir/$objdir/$old_library"
	    old_convenience="$old_convenience $ladir/$objdir/$old_library"
	    tmp_libs=
	    for deplib in $dependency_libs; do
	      deplibs="$deplib $deplibs"
              if test "X$duplicate_deps" = "Xyes" ; then
	        case "$tmp_libs " in
	        *" $deplib "*) specialdeplibs="$specialdeplibs $deplib" ;;
	        esac
              fi
	      tmp_libs="$tmp_libs $deplib"
	    done
	  elif test $linkmode != prog && test $linkmode != lib; then
	    $echo "$modename: \`$lib' is not a convenience library" 1>&2
	    exit 1
	  fi
	  continue
	fi # $pass = conv

	# Get the name of the library we link against.
	linklib=
	for l in $old_library $library_names; do
	  linklib="$l"
	done
	if test -z "$linklib"; then
	  $echo "$modename: cannot find name of link library for \`$lib'" 1>&2
	  exit 1
	fi

	# This library was specified with -dlopen.
	if test $pass = dlopen; then
	  if test -z "$libdir"; then
	    $echo "$modename: cannot -dlopen a convenience library: \`$lib'" 1>&2
	    exit 1
	  fi
	  if test -z "$dlname" || test "$dlopen_support" != yes || test "$build_libtool_libs" = no; then
	    # If there is no dlname, no dlopen support or we're linking
	    # statically, we need to preload.
	    dlprefiles="$dlprefiles $lib"
	  else
	    newdlfiles="$newdlfiles $lib"
	  fi
	  continue
	fi # $pass = dlopen

	# We need an absolute path.
	case $ladir in
	[\\/]* | [A-Za-z]:[\\/]*) abs_ladir="$ladir" ;;
	*)
	  abs_ladir=`cd "$ladir" && pwd`
	  if test -z "$abs_ladir"; then
	    $echo "$modename: warning: cannot determine absolute directory name of \`$ladir'" 1>&2
	    $echo "$modename: passing it literally to the linker, although it might fail" 1>&2
	    abs_ladir="$ladir"
	  fi
	  ;;
	esac
	laname=`$echo "X$lib" | $Xsed -e 's%^.*/%%'`

	# Find the relevant object directory and library name.
	if test "X$installed" = Xyes; then
	  if test ! -f "$libdir/$linklib" && test -f "$abs_ladir/$linklib"; then
	    $echo "$modename: warning: library \`$lib' was moved." 1>&2
	    dir="$ladir"
	    absdir="$abs_ladir"
	    libdir="$abs_ladir"
	  else
	    dir="$libdir"
	    absdir="$libdir"
	  fi
	else
	  dir="$ladir/$objdir"
	  absdir="$abs_ladir/$objdir"
	  # Remove this search path later
	  notinst_path="$notinst_path $abs_ladir"
	fi # $installed = yes
	name=`$echo "X$laname" | $Xsed -e 's/\.la$//' -e 's/^lib//'`

	# This library was specified with -dlpreopen.
	if test $pass = dlpreopen; then
	  if test -z "$libdir"; then
	    $echo "$modename: cannot -dlpreopen a convenience library: \`$lib'" 1>&2
	    exit 1
	  fi
	  # Prefer using a static library (so that no silly _DYNAMIC symbols
	  # are required to link).
	  if test -n "$old_library"; then
	    newdlprefiles="$newdlprefiles $dir/$old_library"
	  # Otherwise, use the dlname, so that lt_dlopen finds it.
	  elif test -n "$dlname"; then
	    newdlprefiles="$newdlprefiles $dir/$dlname"
	  else
	    newdlprefiles="$newdlprefiles $dir/$linklib"
	  fi
	fi # $pass = dlpreopen

	if test -z "$libdir"; then
	  # Link the convenience library
	  if test $linkmode = lib; then
	    deplibs="$dir/$old_library $deplibs"
	  elif test "$linkmode,$pass" = "prog,link"; then
	    compile_deplibs="$dir/$old_library $compile_deplibs"
	    finalize_deplibs="$dir/$old_library $finalize_deplibs"
	  else
	    deplibs="$lib $deplibs"
	  fi
	  continue
	fi

	if test $linkmode = prog && test $pass != link; then
	  newlib_search_path="$newlib_search_path $ladir"
	  deplibs="$lib $deplibs"

	  linkalldeplibs=no
	  if test "$link_all_deplibs" != no || test -z "$library_names" ||
	     test "$build_libtool_libs" = no; then
	    linkalldeplibs=yes
	  fi

	  tmp_libs=
	  for deplib in $dependency_libs; do
	    case $deplib in
	    -L*) newlib_search_path="$newlib_search_path "`$echo "X$deplib" | $Xsed -e 's/^-L//'`;; ### testsuite: skip nested quoting test
	    esac
	    # Need to link against all dependency_libs?
	    if test $linkalldeplibs = yes; then
	      deplibs="$deplib $deplibs"
	    else
	      # Need to hardcode shared library paths
	      # or/and link against static libraries
	      newdependency_libs="$deplib $newdependency_libs"
	    fi
	    if test "X$duplicate_deps" = "Xyes" ; then
	      case "$tmp_libs " in
	      *" $deplib "*) specialdeplibs="$specialdeplibs $deplib" ;;
	      esac
	    fi
	    tmp_libs="$tmp_libs $deplib"
	  done # for deplib
	  continue
	fi # $linkmode = prog...

	link_static=no # Whether the deplib will be linked statically
	if test -n "$library_names" &&
	   { test "$prefer_static_libs" = no || test -z "$old_library"; }; then
	  # Link against this shared library

	  if test "$linkmode,$pass" = "prog,link" ||
	   { test $linkmode = lib && test $hardcode_into_libs = yes; }; then
	    # Hardcode the library path.
	    # Skip directories that are in the system default run-time
	    # search path.
	    case " $sys_lib_dlsearch_path " in
	    *" $absdir "*) ;;
	    *)
	      case "$compile_rpath " in
	      *" $absdir "*) ;;
	      *) compile_rpath="$compile_rpath $absdir"
	      esac
	      ;;
	    esac
	    case " $sys_lib_dlsearch_path " in
	    *" $libdir "*) ;;
	    *)
	      case "$finalize_rpath " in
	      *" $libdir "*) ;;
	      *) finalize_rpath="$finalize_rpath $libdir"
	      esac
	      ;;
	    esac
	    if test $linkmode = prog; then
	      # We need to hardcode the library path
	      if test -n "$shlibpath_var"; then
		# Make sure the rpath contains only unique directories.
		case "$temp_rpath " in
		*" $dir "*) ;;
		*" $absdir "*) ;;
		*) temp_rpath="$temp_rpath $dir" ;;
		esac
	      fi
	    fi
	  fi # $linkmode,$pass = prog,link...

	  if test "$alldeplibs" = yes &&
	     { test "$deplibs_check_method" = pass_all ||
	       { test "$build_libtool_libs" = yes &&
		 test -n "$library_names"; }; }; then
	    # We only need to search for static libraries
	    continue
	  fi

	  if test "$installed" = no; then
	    notinst_deplibs="$notinst_deplibs $lib"
	    need_relink=yes
	  fi

	  if test -n "$old_archive_from_expsyms_cmds"; then
	    # figure out the soname
	    set dummy $library_names
	    realname="$2"
	    shift; shift
	    libname=`eval \\$echo \"$libname_spec\"`
	    # use dlname if we got it. it's perfectly good, no?
	    if test -n "$dlname"; then
	      soname="$dlname"
	    elif test -n "$soname_spec"; then
	      # bleh windows
	      case $host in
	      *cygwin*)
		major=`expr $current - $age`
		versuffix="-$major"
		;;
	      esac
	      eval soname=\"$soname_spec\"
	    else
	      soname="$realname"
	    fi

	    # Make a new name for the extract_expsyms_cmds to use
	    soroot="$soname"
	    soname=`echo $soroot | ${SED} -e 's/^.*\///'`
	    newlib="libimp-`echo $soname | ${SED} 's/^lib//;s/\.dll$//'`.a"

	    # If the library has no export list, then create one now
	    if test -f "$output_objdir/$soname-def"; then :
	    else
	      $show "extracting exported symbol list from \`$soname'"
	      save_ifs="$IFS"; IFS='~'
	      eval cmds=\"$extract_expsyms_cmds\"
	      for cmd in $cmds; do
		IFS="$save_ifs"
		$show "$cmd"
		$run eval "$cmd" || exit $?
	      done
	      IFS="$save_ifs"
	    fi

	    # Create $newlib
	    if test -f "$output_objdir/$newlib"; then :; else
	      $show "generating import library for \`$soname'"
	      save_ifs="$IFS"; IFS='~'
	      eval cmds=\"$old_archive_from_expsyms_cmds\"
	      for cmd in $cmds; do
		IFS="$save_ifs"
		$show "$cmd"
		$run eval "$cmd" || exit $?
	      done
	      IFS="$save_ifs"
	    fi
	    # make sure the library variables are pointing to the new library
	    dir=$output_objdir
	    linklib=$newlib
	  fi # test -n $old_archive_from_expsyms_cmds

	  if test $linkmode = prog || test "$mode" != relink; then
	    add_shlibpath=
	    add_dir=
	    add=
	    lib_linked=yes
	    case $hardcode_action in
	    immediate | unsupported)
	      if test "$hardcode_direct" = no; then
		add="$dir/$linklib"
	      elif test "$hardcode_minus_L" = no; then
		case $host in
		*-*-sunos*) add_shlibpath="$dir" ;;
		esac
		add_dir="-L$dir"
		add="-l$name"
	      elif test "$hardcode_shlibpath_var" = no; then
		add_shlibpath="$dir"
		add="-l$name"
	      else
		lib_linked=no
	      fi
	      ;;
	    relink)
	      if test "$hardcode_direct" = yes; then
		add="$dir/$linklib"
	      elif test "$hardcode_minus_L" = yes; then
		add_dir="-L$dir"
		add="-l$name"
	      elif test "$hardcode_shlibpath_var" = yes; then
		add_shlibpath="$dir"
		add="-l$name"
	      else
		lib_linked=no
	      fi
	      ;;
	    *) lib_linked=no ;;
	    esac

	    if test "$lib_linked" != yes; then
	      $echo "$modename: configuration error: unsupported hardcode properties"
	      exit 1
	    fi

	    if test -n "$add_shlibpath"; then
	      case :$compile_shlibpath: in
	      *":$add_shlibpath:"*) ;;
	      *) compile_shlibpath="$compile_shlibpath$add_shlibpath:" ;;
	      esac
	    fi
	    if test $linkmode = prog; then
	      test -n "$add_dir" && compile_deplibs="$add_dir $compile_deplibs"
	      test -n "$add" && compile_deplibs="$add $compile_deplibs"
	    else
	      test -n "$add_dir" && deplibs="$add_dir $deplibs"
	      test -n "$add" && deplibs="$add $deplibs"
	      if test "$hardcode_direct" != yes && \
		 test "$hardcode_minus_L" != yes && \
		 test "$hardcode_shlibpath_var" = yes; then
		case :$finalize_shlibpath: in
		*":$libdir:"*) ;;
		*) finalize_shlibpath="$finalize_shlibpath$libdir:" ;;
		esac
	      fi
	    fi
	  fi

	  if test $linkmode = prog || test "$mode" = relink; then
	    add_shlibpath=
	    add_dir=
	    add=
	    # Finalize command for both is simple: just hardcode it.
	    if test "$hardcode_direct" = yes; then
	      add="$libdir/$linklib"
	    elif test "$hardcode_minus_L" = yes; then
	      add_dir="-L$libdir"
	      add="-l$name"
	    elif test "$hardcode_shlibpath_var" = yes; then
	      case :$finalize_shlibpath: in
	      *":$libdir:"*) ;;
	      *) finalize_shlibpath="$finalize_shlibpath$libdir:" ;;
	      esac
	      add="-l$name"
	    else
	      # We cannot seem to hardcode it, guess we'll fake it.
	      if test "X$installed" = Xyes; then
	        add_dir="-L$libdir"
	      else
	        add_dir="-L$DESTDIR$libdir"
	      fi
	      add="-l$name"
	    fi

	    if test $linkmode = prog; then
	      test -n "$add_dir" && finalize_deplibs="$add_dir $finalize_deplibs"
	      test -n "$add" && finalize_deplibs="$add $finalize_deplibs"
	    else
	      test -n "$add_dir" && deplibs="$add_dir $deplibs"
	      test -n "$add" && deplibs="$add $deplibs"
	    fi
	  fi
	elif test $linkmode = prog; then
	  if test "$alldeplibs" = yes &&
	     { test "$deplibs_check_method" = pass_all ||
	       { test "$build_libtool_libs" = yes &&
		 test -n "$library_names"; }; }; then
	    # We only need to search for static libraries
	    continue
	  fi

	  # Try to link the static library
	  # Here we assume that one of hardcode_direct or hardcode_minus_L
	  # is not unsupported.  This is valid on all known static and
	  # shared platforms.
	  if test "$hardcode_direct" != unsupported; then
	    test -n "$old_library" && linklib="$old_library"
	    compile_deplibs="$dir/$linklib $compile_deplibs"
	    finalize_deplibs="$dir/$linklib $finalize_deplibs"
	  else
	    compile_deplibs="-l$name -L$dir $compile_deplibs"
	    finalize_deplibs="-l$name -L$dir $finalize_deplibs"
	  fi
	elif test "$build_libtool_libs" = yes; then
	  # Not a shared library
	  if test "$deplibs_check_method" != pass_all; then
	    # We're trying link a shared library against a static one
	    # but the system doesn't support it.

	    # Just print a warning and add the library to dependency_libs so
	    # that the program can be linked against the static library.
	    echo
	    echo "*** Warning: This system can not link to static lib archive $lib."
	    echo "*** I have the capability to make that library automatically link in when"
	    echo "*** you link to this library.  But I can only do this if you have a"
	    echo "*** shared version of the library, which you do not appear to have."
	    if test "$module" = yes; then
	      echo "*** But as you try to build a module library, libtool will still create "
	      echo "*** a static module, that should work as long as the dlopening application"
	      echo "*** is linked with the -dlopen flag to resolve symbols at runtime."
	      if test -z "$global_symbol_pipe"; then
		echo
		echo "*** However, this would only work if libtool was able to extract symbol"
		echo "*** lists from a program, using \`nm' or equivalent, but libtool could"
		echo "*** not find such a program.  So, this module is probably useless."
		echo "*** \`nm' from GNU binutils and a full rebuild may help."
	      fi
	      if test "$build_old_libs" = no; then
		build_libtool_libs=module
		build_old_libs=yes
	      else
		build_libtool_libs=no
	      fi
	    fi
	  else
	    convenience="$convenience $dir/$old_library"
	    old_convenience="$old_convenience $dir/$old_library"
	    deplibs="$dir/$old_library $deplibs"
	    link_static=yes
	  fi
	fi # link shared/static library?

	if test $linkmode = lib; then
	  if test -n "$dependency_libs" &&
	     { test $hardcode_into_libs != yes || test $build_old_libs = yes ||
	       test $link_static = yes; }; then
	    # Extract -R from dependency_libs
	    temp_deplibs=
	    for libdir in $dependency_libs; do
	      case $libdir in
	      -R*) temp_xrpath=`$echo "X$libdir" | $Xsed -e 's/^-R//'`
		   case " $xrpath " in
		   *" $temp_xrpath "*) ;;
		   *) xrpath="$xrpath $temp_xrpath";;
		   esac;;
	      *) temp_deplibs="$temp_deplibs $libdir";;
	      esac
	    done
	    dependency_libs="$temp_deplibs"
	  fi

	  newlib_search_path="$newlib_search_path $absdir"
	  # Link against this library
	  test "$link_static" = no && newdependency_libs="$abs_ladir/$laname $newdependency_libs"
	  # ... and its dependency_libs
	  tmp_libs=
	  for deplib in $dependency_libs; do
	    newdependency_libs="$deplib $newdependency_libs"
	    if test "X$duplicate_deps" = "Xyes" ; then
	      case "$tmp_libs " in
	      *" $deplib "*) specialdeplibs="$specialdeplibs $deplib" ;;
	      esac
	    fi
	    tmp_libs="$tmp_libs $deplib"
	  done

	  if test $link_all_deplibs != no; then
	    # Add the search paths of all dependency libraries
	    for deplib in $dependency_libs; do
	      case $deplib in
	      -L*) path="$deplib" ;;
	      *.la)
		dir=`$echo "X$deplib" | $Xsed -e 's%/[^/]*$%%'`
		test "X$dir" = "X$deplib" && dir="."
		# We need an absolute path.
		case $dir in
		[\\/]* | [A-Za-z]:[\\/]*) absdir="$dir" ;;
		*)
		  absdir=`cd "$dir" && pwd`
		  if test -z "$absdir"; then
		    $echo "$modename: warning: cannot determine absolute directory name of \`$dir'" 1>&2
		    absdir="$dir"
		  fi
		  ;;
		esac
		if grep "^installed=no" $deplib > /dev/null; then
		  path="-L$absdir/$objdir"
		else
		  eval libdir=`${SED} -n -e 's/^libdir=\(.*\)$/\1/p' $deplib`
		  if test -z "$libdir"; then
		    $echo "$modename: \`$deplib' is not a valid libtool archive" 1>&2
		    exit 1
		  fi
		  if test "$absdir" != "$libdir"; then
		    $echo "$modename: warning: \`$deplib' seems to be moved" 1>&2
		  fi
		  path="-L$absdir"
		fi
		;;
	      *) continue ;;
	      esac
	      case " $deplibs " in
	      *" $path "*) ;;
	      *) deplibs="$deplibs $path" ;;
	      esac
	    done
	  fi # link_all_deplibs != no
	fi # linkmode = lib
      done # for deplib in $libs
      if test $pass = dlpreopen; then
	# Link the dlpreopened libraries before other libraries
	for deplib in $save_deplibs; do
	  deplibs="$deplib $deplibs"
	done
      fi
      if test $pass != dlopen; then
	test $pass != scan && dependency_libs="$newdependency_libs"
	if test $pass != conv; then
	  # Make sure lib_search_path contains only unique directories.
	  lib_search_path=
	  for dir in $newlib_search_path; do
	    case "$lib_search_path " in
	    *" $dir "*) ;;
	    *) lib_search_path="$lib_search_path $dir" ;;
	    esac
	  done
	  newlib_search_path=
	fi

	if test "$linkmode,$pass" != "prog,link"; then
	  vars="deplibs"
	else
	  vars="compile_deplibs finalize_deplibs"
	fi
	for var in $vars dependency_libs; do
	  # Add libraries to $var in reverse order
	  eval tmp_libs=\"\$$var\"
	  new_libs=
	  for deplib in $tmp_libs; do
	    case $deplib in
	    -L*) new_libs="$deplib $new_libs" ;;
	    *)
	      case " $specialdeplibs " in
	      *" $deplib "*) new_libs="$deplib $new_libs" ;;
	      *)
		case " $new_libs " in
		*" $deplib "*) ;;
		*) new_libs="$deplib $new_libs" ;;
		esac
		;;
	      esac
	      ;;
	    esac
	  done
	  tmp_libs=
	  for deplib in $new_libs; do
	    case $deplib in
	    -L*)
	      case " $tmp_libs " in
	      *" $deplib "*) ;;
	      *) tmp_libs="$tmp_libs $deplib" ;;
	      esac
	      ;;
	    *) tmp_libs="$tmp_libs $deplib" ;;
	    esac
	  done
	  eval $var=\"$tmp_libs\"
	done # for var
      fi
      if test "$pass" = "conv" &&
       { test "$linkmode" = "lib" || test "$linkmode" = "prog"; }; then
	libs="$deplibs" # reset libs
	deplibs=
      fi
    done # for pass
    if test $linkmode = prog; then
      dlfiles="$newdlfiles"
      dlprefiles="$newdlprefiles"
    fi

    case $linkmode in
    oldlib)
      if test -n "$dlfiles$dlprefiles" || test "$dlself" != no; then
	$echo "$modename: warning: \`-dlopen' is ignored for archives" 1>&2
      fi

      if test -n "$rpath"; then
	$echo "$modename: warning: \`-rpath' is ignored for archives" 1>&2
      fi

      if test -n "$xrpath"; then
	$echo "$modename: warning: \`-R' is ignored for archives" 1>&2
      fi

      if test -n "$vinfo"; then
	$echo "$modename: warning: \`-version-info' is ignored for archives" 1>&2
      fi

      if test -n "$release"; then
	$echo "$modename: warning: \`-release' is ignored for archives" 1>&2
      fi

      if test -n "$export_symbols" || test -n "$export_symbols_regex"; then
	$echo "$modename: warning: \`-export-symbols' is ignored for archives" 1>&2
      fi

      # Now set the variables for building old libraries.
      build_libtool_libs=no
      oldlibs="$output"
      objs="$objs$old_deplibs"
      ;;

    lib)
      # Make sure we only generate libraries of the form `libNAME.la'.
      case $outputname in
      lib*)
	name=`$echo "X$outputname" | $Xsed -e 's/\.la$//' -e 's/^lib//'`
	eval libname=\"$libname_spec\"
	;;
      *)
	if test "$module" = no; then
	  $echo "$modename: libtool library \`$output' must begin with \`lib'" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	fi
	if test "$need_lib_prefix" != no; then
	  # Add the "lib" prefix for modules if required
	  name=`$echo "X$outputname" | $Xsed -e 's/\.la$//'`
	  eval libname=\"$libname_spec\"
	else
	  libname=`$echo "X$outputname" | $Xsed -e 's/\.la$//'`
	fi
	;;
      esac

      if test -n "$objs"; then
	if test "$deplibs_check_method" != pass_all; then
	  $echo "$modename: cannot build libtool library \`$output' from non-libtool objects on this host:$objs" 2>&1
	  exit 1
	else
	  echo
	  echo "*** Warning: Linking the shared library $output against the non-libtool"
	  echo "*** objects $objs is not portable!"
	  libobjs="$libobjs $objs"
	fi
      fi

      if test "$dlself" != no; then
	$echo "$modename: warning: \`-dlopen self' is ignored for libtool libraries" 1>&2
      fi

      set dummy $rpath
      if test $# -gt 2; then
	$echo "$modename: warning: ignoring multiple \`-rpath's for a libtool library" 1>&2
      fi
      install_libdir="$2"

      oldlibs=
      if test -z "$rpath"; then
	if test "$build_libtool_libs" = yes; then
	  # Building a libtool convenience library.
	  libext=al
	  oldlibs="$output_objdir/$libname.$libext $oldlibs"
	  build_libtool_libs=convenience
	  build_old_libs=yes
	fi

	if test -n "$vinfo"; then
	  $echo "$modename: warning: \`-version-info' is ignored for convenience libraries" 1>&2
	fi

	if test -n "$release"; then
	  $echo "$modename: warning: \`-release' is ignored for convenience libraries" 1>&2
	fi
      else

	# Parse the version information argument.
	save_ifs="$IFS"; IFS=':'
	set dummy $vinfo 0 0 0
	IFS="$save_ifs"

	if test -n "$8"; then
	  $echo "$modename: too many parameters to \`-version-info'" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	fi

	current="$2"
	revision="$3"
	age="$4"

	# Check that each of the things are valid numbers.
	case $current in
	0 | [1-9] | [1-9][0-9] | [1-9][0-9][0-9]) ;;
	*)
	  $echo "$modename: CURRENT \`$current' is not a nonnegative integer" 1>&2
	  $echo "$modename: \`$vinfo' is not valid version information" 1>&2
	  exit 1
	  ;;
	esac

	case $revision in
	0 | [1-9] | [1-9][0-9] | [1-9][0-9][0-9]) ;;
	*)
	  $echo "$modename: REVISION \`$revision' is not a nonnegative integer" 1>&2
	  $echo "$modename: \`$vinfo' is not valid version information" 1>&2
	  exit 1
	  ;;
	esac

	case $age in
	0 | [1-9] | [1-9][0-9] | [1-9][0-9][0-9]) ;;
	*)
	  $echo "$modename: AGE \`$age' is not a nonnegative integer" 1>&2
	  $echo "$modename: \`$vinfo' is not valid version information" 1>&2
	  exit 1
	  ;;
	esac

	if test $age -gt $current; then
	  $echo "$modename: AGE \`$age' is greater than the current interface number \`$current'" 1>&2
	  $echo "$modename: \`$vinfo' is not valid version information" 1>&2
	  exit 1
	fi

	# Calculate the version variables.
	major=
	versuffix=
	verstring=
	case $version_type in
	none) ;;

	darwin)
	  # Like Linux, but with the current version available in
	  # verstring for coding it into the library header
	  major=.`expr $current - $age`
	  versuffix="$major.$age.$revision"
	  # Darwin ld doesn't like 0 for these options...
	  minor_current=`expr $current + 1`
	  verstring="-compatibility_version $minor_current -current_version $minor_current.$revision"
	  ;;

	freebsd-aout)
	  major=".$current"
	  versuffix=".$current.$revision";
	  ;;

	freebsd-elf)
	  major=".$current"
	  versuffix=".$current";
	  ;;

	irix | nonstopux)
	  major=`expr $current - $age + 1`

	  case $version_type in
	    nonstopux) verstring_prefix=nonstopux ;;
	    *)         verstring_prefix=sgi ;;
	  esac
	  verstring="$verstring_prefix$major.$revision"

	  # Add in all the interfaces that we are compatible with.
	  loop=$revision
	  while test $loop != 0; do
	    iface=`expr $revision - $loop`
	    loop=`expr $loop - 1`
	    verstring="$verstring_prefix$major.$iface:$verstring"
	  done

	  # Before this point, $major must not contain `.'.
	  major=.$major
	  versuffix="$major.$revision"
	  ;;

	linux)
	  major=.`expr $current - $age`
	  versuffix="$major.$age.$revision"
	  ;;

	osf)
	  major=.`expr $current - $age`
	  versuffix=".$current.$age.$revision"
	  verstring="$current.$age.$revision"

	  # Add in all the interfaces that we are compatible with.
	  loop=$age
	  while test $loop != 0; do
	    iface=`expr $current - $loop`
	    loop=`expr $loop - 1`
	    verstring="$verstring:${iface}.0"
	  done

	  # Make executables depend on our current version.
	  verstring="$verstring:${current}.0"
	  ;;

	sunos)
	  major=".$current"
	  versuffix=".$current.$revision"
	  ;;

	windows)
	  # Use '-' rather than '.', since we only want one
	  # extension on DOS 8.3 filesystems.
	  major=`expr $current - $age`
	  versuffix="-$major"
	  ;;

	*)
	  $echo "$modename: unknown library version type \`$version_type'" 1>&2
	  echo "Fatal configuration error.  See the $PACKAGE docs for more information." 1>&2
	  exit 1
	  ;;
	esac

	# Clear the version info if we defaulted, and they specified a release.
	if test -z "$vinfo" && test -n "$release"; then
	  major=
	  verstring="0.0"
	  case $version_type in
	  darwin)
	    # we can't check for "0.0" in archive_cmds due to quoting
	    # problems, so we reset it completely
	    verstring=""
	    ;;
	  *)
	    verstring="0.0"
	    ;;
	  esac
	  if test "$need_version" = no; then
	    versuffix=
	  else
	    versuffix=".0.0"
	  fi
	fi

	# Remove version info from name if versioning should be avoided
	if test "$avoid_version" = yes && test "$need_version" = no; then
	  major=
	  versuffix=
	  verstring=""
	fi

	# Check to see if the archive will have undefined symbols.
	if test "$allow_undefined" = yes; then
	  if test "$allow_undefined_flag" = unsupported; then
	    $echo "$modename: warning: undefined symbols not allowed in $host shared libraries" 1>&2
	    build_libtool_libs=no
	    build_old_libs=yes
	  fi
	else
	  # Don't allow undefined symbols.
	  allow_undefined_flag="$no_undefined_flag"
	fi
      fi

      if test "$mode" != relink; then
	# Remove our outputs.
	$show "${rm}r $output_objdir/$outputname $output_objdir/$libname.* $output_objdir/${libname}${release}.*"
	$run ${rm}r $output_objdir/$outputname $output_objdir/$libname.* $output_objdir/${libname}${release}.*
      fi

      # Now set the variables for building old libraries.
      if test "$build_old_libs" = yes && test "$build_libtool_libs" != convenience ; then
	oldlibs="$oldlibs $output_objdir/$libname.$libext"

	# Transform .lo files to .o files.
	oldobjs="$objs "`$echo "X$libobjs" | $SP2NL | $Xsed -e '/\.'${libext}'$/d' -e "$lo2o" | $NL2SP`
      fi

      # Eliminate all temporary directories.
      for path in $notinst_path; do
	lib_search_path=`echo "$lib_search_path " | ${SED} -e 's% $path % %g'`
	deplibs=`echo "$deplibs " | ${SED} -e 's% -L$path % %g'`
	dependency_libs=`echo "$dependency_libs " | ${SED} -e 's% -L$path % %g'`
      done

      if test -n "$xrpath"; then
	# If the user specified any rpath flags, then add them.
	temp_xrpath=
	for libdir in $xrpath; do
	  temp_xrpath="$temp_xrpath -R$libdir"
	  case "$finalize_rpath " in
	  *" $libdir "*) ;;
	  *) finalize_rpath="$finalize_rpath $libdir" ;;
	  esac
	done
	if test $hardcode_into_libs != yes || test $build_old_libs = yes; then
	  dependency_libs="$temp_xrpath $dependency_libs"
	fi
      fi

      # Make sure dlfiles contains only unique files that won't be dlpreopened
      old_dlfiles="$dlfiles"
      dlfiles=
      for lib in $old_dlfiles; do
	case " $dlprefiles $dlfiles " in
	*" $lib "*) ;;
	*) dlfiles="$dlfiles $lib" ;;
	esac
      done

      # Make sure dlprefiles contains only unique files
      old_dlprefiles="$dlprefiles"
      dlprefiles=
      for lib in $old_dlprefiles; do
	case "$dlprefiles " in
	*" $lib "*) ;;
	*) dlprefiles="$dlprefiles $lib" ;;
	esac
      done

      if test "$build_libtool_libs" = yes; then
	if test -n "$rpath"; then
	  case $host in
	  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2* | *-*-beos*)
	    # these systems don't actually have a c library (as such)!
	    ;;
	  *-*-rhapsody* | *-*-darwin1.[012])
	    # Rhapsody C library is in the System framework
	    deplibs="$deplibs -framework System"
	    ;;
	  *-*-netbsd*)
	    # Don't link with libc until the a.out ld.so is fixed.
	    ;;
	  *-*-openbsd* | *-*-freebsd*)
	    # Do not include libc due to us having libc/libc_r.
	    ;;
	  *)
	    # Add libc to deplibs on all other systems if necessary.
	    if test $build_libtool_need_lc = "yes"; then
	      deplibs="$deplibs -lc"
	    fi
	    ;;
	  esac
	fi

	# Transform deplibs into only deplibs that can be linked in shared.
	name_save=$name
	libname_save=$libname
	release_save=$release
	versuffix_save=$versuffix
	major_save=$major
	# I'm not sure if I'm treating the release correctly.  I think
	# release should show up in the -l (ie -lgmp5) so we don't want to
	# add it in twice.  Is that correct?
	release=""
	versuffix=""
	major=""
	newdeplibs=
	droppeddeps=no
	case $deplibs_check_method in
	pass_all)
	  # Don't check for shared/static.  Everything works.
	  # This might be a little naive.  We might want to check
	  # whether the library exists or not.  But this is on
	  # osf3 & osf4 and I'm not really sure... Just
	  # implementing what was already the behaviour.
	  newdeplibs=$deplibs
	  ;;
	test_compile)
	  # This code stresses the "libraries are programs" paradigm to its
	  # limits. Maybe even breaks it.  We compile a program, linking it
	  # against the deplibs as a proxy for the library.  Then we can check
	  # whether they linked in statically or dynamically with ldd.
	  $rm conftest.c
	  cat > conftest.c <<EOF
	  int main() { return 0; }
EOF
	  $rm conftest
	  $CC -o conftest conftest.c $deplibs
	  if test $? -eq 0 ; then
	    ldd_output=`ldd conftest`
	    for i in $deplibs; do
	      name="`expr $i : '-l\(.*\)'`"
	      # If $name is empty we are operating on a -L argument.
	      if test -n "$name" && test "$name" != "0"; then
		libname=`eval \\$echo \"$libname_spec\"`
		deplib_matches=`eval \\$echo \"$library_names_spec\"`
		set dummy $deplib_matches
		deplib_match=$2
		if test `expr "$ldd_output" : ".*$deplib_match"` -ne 0 ; then
		  newdeplibs="$newdeplibs $i"
		else
		  droppeddeps=yes
		  echo
		  echo "*** Warning: dynamic linker does not accept needed library $i."
		  echo "*** I have the capability to make that library automatically link in when"
		  echo "*** you link to this library.  But I can only do this if you have a"
		  echo "*** shared version of the library, which I believe you do not have"
		  echo "*** because a test_compile did reveal that the linker did not use it for"
		  echo "*** its dynamic dependency list that programs get resolved with at runtime."
		fi
	      else
		newdeplibs="$newdeplibs $i"
	      fi
	    done
	  else
	    # Error occured in the first compile.  Let's try to salvage
	    # the situation: Compile a separate program for each library.
	    for i in $deplibs; do
	      name="`expr $i : '-l\(.*\)'`"
	     # If $name is empty we are operating on a -L argument.
	      if test -n "$name" && test "$name" != "0"; then
		$rm conftest
		$CC -o conftest conftest.c $i
		# Did it work?
		if test $? -eq 0 ; then
		  ldd_output=`ldd conftest`
		  libname=`eval \\$echo \"$libname_spec\"`
		  deplib_matches=`eval \\$echo \"$library_names_spec\"`
		  set dummy $deplib_matches
		  deplib_match=$2
		  if test `expr "$ldd_output" : ".*$deplib_match"` -ne 0 ; then
		    newdeplibs="$newdeplibs $i"
		  else
		    droppeddeps=yes
		    echo
		    echo "*** Warning: dynamic linker does not accept needed library $i."
		    echo "*** I have the capability to make that library automatically link in when"
		    echo "*** you link to this library.  But I can only do this if you have a"
		    echo "*** shared version of the library, which you do not appear to have"
		    echo "*** because a test_compile did reveal that the linker did not use this one"
		    echo "*** as a dynamic dependency that programs can get resolved with at runtime."
		  fi
		else
		  droppeddeps=yes
		  echo
		  echo "*** Warning!  Library $i is needed by this library but I was not able to"
		  echo "***  make it link in!  You will probably need to install it or some"
		  echo "*** library that it depends on before this library will be fully"
		  echo "*** functional.  Installing it before continuing would be even better."
		fi
	      else
		newdeplibs="$newdeplibs $i"
	      fi
	    done
	  fi
	  ;;
	file_magic*)
	  set dummy $deplibs_check_method
	  file_magic_regex=`expr "$deplibs_check_method" : "$2 \(.*\)"`
	  for a_deplib in $deplibs; do
	    name="`expr $a_deplib : '-l\(.*\)'`"
	    # If $name is empty we are operating on a -L argument.
	    if test -n "$name" && test "$name" != "0"; then
	      libname=`eval \\$echo \"$libname_spec\"`
	      for i in $lib_search_path $sys_lib_search_path $shlib_search_path; do
		    potential_libs=`ls $i/$libname[.-]* 2>/dev/null`
		    for potent_lib in $potential_libs; do
		      # Follow soft links.
		      if ls -lLd "$potent_lib" 2>/dev/null \
			 | grep " -> " >/dev/null; then
			continue
		      fi
		      # The statement above tries to avoid entering an
		      # endless loop below, in case of cyclic links.
		      # We might still enter an endless loop, since a link
		      # loop can be closed while we follow links,
		      # but so what?
		      potlib="$potent_lib"
		      while test -h "$potlib" 2>/dev/null; do
			potliblink=`ls -ld $potlib | ${SED} 's/.* -> //'`
			case $potliblink in
			[\\/]* | [A-Za-z]:[\\/]*) potlib="$potliblink";;
			*) potlib=`$echo "X$potlib" | $Xsed -e 's,[^/]*$,,'`"$potliblink";;
			esac
		      done
		      # It is ok to link against an archive when
		      # building a shared library.
		      if $AR -t $potlib > /dev/null 2>&1; then
		        newdeplibs="$newdeplibs $a_deplib"
		        a_deplib=""
		        break 2
		      fi
		      if eval $file_magic_cmd \"\$potlib\" 2>/dev/null \
			 | ${SED} 10q \
			 | egrep "$file_magic_regex" > /dev/null; then
			newdeplibs="$newdeplibs $a_deplib"
			a_deplib=""
			break 2
		      fi
		    done
	      done
	      if test -n "$a_deplib" ; then
		droppeddeps=yes
		echo
		echo "*** Warning: linker path does not have real file for library $a_deplib."
		echo "*** I have the capability to make that library automatically link in when"
		echo "*** you link to this library.  But I can only do this if you have a"
		echo "*** shared version of the library, which you do not appear to have"
		echo "*** because I did check the linker path looking for a file starting"
		if test -z "$potlib" ; then
		  echo "*** with $libname but no candidates were found. (...for file magic test)"
		else
		  echo "*** with $libname and none of the candidates passed a file format test"
		  echo "*** using a file magic. Last file checked: $potlib"
		fi
	      fi
	    else
	      # Add a -L argument.
	      newdeplibs="$newdeplibs $a_deplib"
	    fi
	  done # Gone through all deplibs.
	  ;;
	match_pattern*)
	  set dummy $deplibs_check_method
	  match_pattern_regex=`expr "$deplibs_check_method" : "$2 \(.*\)"`
	  for a_deplib in $deplibs; do
	    name="`expr $a_deplib : '-l\(.*\)'`"
	    # If $name is empty we are operating on a -L argument.
	    if test -n "$name" && test "$name" != "0"; then
	      libname=`eval \\$echo \"$libname_spec\"`
	      for i in $lib_search_path $sys_lib_search_path $shlib_search_path; do
		potential_libs=`ls $i/$libname[.-]* 2>/dev/null`
		for potent_lib in $potential_libs; do
		  potlib="$potent_lib" # see symlink-check below in file_magic test
		  if eval echo \"$potent_lib\" 2>/dev/null \
		      | ${SED} 10q \
		      | egrep "$match_pattern_regex" > /dev/null; then
		    newdeplibs="$newdeplibs $a_deplib"
		    a_deplib=""
		    break 2
		  fi
		done
	      done
	      if test -n "$a_deplib" ; then
		droppeddeps=yes
		echo
		echo "*** Warning: linker path does not have real file for library $a_deplib."
		echo "*** I have the capability to make that library automatically link in when"
		echo "*** you link to this library.  But I can only do this if you have a"
		echo "*** shared version of the library, which you do not appear to have"
		echo "*** because I did check the linker path looking for a file starting"
		if test -z "$potlib" ; then
		  echo "*** with $libname but no candidates were found. (...for regex pattern test)"
		else
		  echo "*** with $libname and none of the candidates passed a file format test"
		  echo "*** using a regex pattern. Last file checked: $potlib"
		fi
	      fi
	    else
	      # Add a -L argument.
	      newdeplibs="$newdeplibs $a_deplib"
	    fi
	  done # Gone through all deplibs.
	  ;;
	none | unknown | *)
	  newdeplibs=""
	  if $echo "X $deplibs" | $Xsed -e 's/ -lc$//' \
	       -e 's/ -[LR][^ ]*//g' -e 's/[ 	]//g' |
	     grep . >/dev/null; then
	    echo
	    if test "X$deplibs_check_method" = "Xnone"; then
	      echo "*** Warning: inter-library dependencies are not supported in this platform."
	    else
	      echo "*** Warning: inter-library dependencies are not known to be supported."
	    fi
	    echo "*** All declared inter-library dependencies are being dropped."
	    droppeddeps=yes
	  fi
	  ;;
	esac
	versuffix=$versuffix_save
	major=$major_save
	release=$release_save
	libname=$libname_save
	name=$name_save

	case $host in
	*-*-rhapsody* | *-*-darwin1.[012])
	  # On Rhapsody replace the C library is the System framework
	  newdeplibs=`$echo "X $newdeplibs" | $Xsed -e 's/ -lc / -framework System /'`
	  ;;
	esac

	if test "$droppeddeps" = yes; then
	  if test "$module" = yes; then
	    echo
	    echo "*** Warning: libtool could not satisfy all declared inter-library"
	    echo "*** dependencies of module $libname.  Therefore, libtool will create"
	    echo "*** a static module, that should work as long as the dlopening"
	    echo "*** application is linked with the -dlopen flag."
	    if test -z "$global_symbol_pipe"; then
	      echo
	      echo "*** However, this would only work if libtool was able to extract symbol"
	      echo "*** lists from a program, using \`nm' or equivalent, but libtool could"
	      echo "*** not find such a program.  So, this module is probably useless."
	      echo "*** \`nm' from GNU binutils and a full rebuild may help."
	    fi
	    if test "$build_old_libs" = no; then
	      oldlibs="$output_objdir/$libname.$libext"
	      build_libtool_libs=module
	      build_old_libs=yes
	    else
	      build_libtool_libs=no
	    fi
	  else
	    echo "*** The inter-library dependencies that have been dropped here will be"
	    echo "*** automatically added whenever a program is linked with this library"
	    echo "*** or is declared to -dlopen it."

	    if test $allow_undefined = no; then
	      echo
	      echo "*** Since this library must not contain undefined symbols,"
	      echo "*** because either the platform does not support them or"
	      echo "*** it was explicitly requested with -no-undefined,"
	      echo "*** libtool will only create a static version of it."
	      if test "$build_old_libs" = no; then
		oldlibs="$output_objdir/$libname.$libext"
		build_libtool_libs=module
		build_old_libs=yes
	      else
		build_libtool_libs=no
	      fi
	    fi
	  fi
	fi
	# Done checking deplibs!
	deplibs=$newdeplibs
      fi

      # All the library-specific variables (install_libdir is set above).
      library_names=
      old_library=
      dlname=

      # Test again, we may have decided not to build it any more
      if test "$build_libtool_libs" = yes; then
	if test $hardcode_into_libs = yes; then
	  # Hardcode the library paths
	  hardcode_libdirs=
	  dep_rpath=
	  rpath="$finalize_rpath"
	  test "$mode" != relink && rpath="$compile_rpath$rpath"
	  for libdir in $rpath; do
	    if test -n "$hardcode_libdir_flag_spec"; then
	      if test -n "$hardcode_libdir_separator"; then
		if test -z "$hardcode_libdirs"; then
		  hardcode_libdirs="$libdir"
		else
		  # Just accumulate the unique libdirs.
		  case $hardcode_libdir_separator$hardcode_libdirs$hardcode_libdir_separator in
		  *"$hardcode_libdir_separator$libdir$hardcode_libdir_separator"*)
		    ;;
		  *)
		    hardcode_libdirs="$hardcode_libdirs$hardcode_libdir_separator$libdir"
		    ;;
		  esac
		fi
	      else
		eval flag=\"$hardcode_libdir_flag_spec\"
		dep_rpath="$dep_rpath $flag"
	      fi
	    elif test -n "$runpath_var"; then
	      case "$perm_rpath " in
	      *" $libdir "*) ;;
	      *) perm_rpath="$perm_rpath $libdir" ;;
	      esac
	    fi
	  done
	  # Substitute the hardcoded libdirs into the rpath.
	  if test -n "$hardcode_libdir_separator" &&
	     test -n "$hardcode_libdirs"; then
	    libdir="$hardcode_libdirs"
	    eval dep_rpath=\"$hardcode_libdir_flag_spec\"
	  fi
	  if test -n "$runpath_var" && test -n "$perm_rpath"; then
	    # We should set the runpath_var.
	    rpath=
	    for dir in $perm_rpath; do
	      rpath="$rpath$dir:"
	    done
	    eval "$runpath_var='$rpath\$$runpath_var'; export $runpath_var"
	  fi
	  test -n "$dep_rpath" && deplibs="$dep_rpath $deplibs"
	fi

	shlibpath="$finalize_shlibpath"
	test "$mode" != relink && shlibpath="$compile_shlibpath$shlibpath"
	if test -n "$shlibpath"; then
	  eval "$shlibpath_var='$shlibpath\$$shlibpath_var'; export $shlibpath_var"
	fi

	# Get the real and link names of the library.
	eval library_names=\"$library_names_spec\"
	set dummy $library_names
	realname="$2"
	shift; shift

	if test -n "$soname_spec"; then
	  eval soname=\"$soname_spec\"
	else
	  soname="$realname"
	fi
	test -z "$dlname" && dlname=$soname

	lib="$output_objdir/$realname"
	for link
	do
	  linknames="$linknames $link"
	done

	# Ensure that we have .o objects for linkers which dislike .lo
	# (e.g. aix) in case we are running --disable-static
	for obj in $libobjs; do
	  xdir=`$echo "X$obj" | $Xsed -e 's%/[^/]*$%%'`
	  if test "X$xdir" = "X$obj"; then
	    xdir="."
	  else
	    xdir="$xdir"
	  fi
	  baseobj=`$echo "X$obj" | $Xsed -e 's%^.*/%%'`
	  oldobj=`$echo "X$baseobj" | $Xsed -e "$lo2o"`
	  if test ! -f $xdir/$oldobj; then
	    $show "(cd $xdir && ${LN_S} $baseobj $oldobj)"
	    $run eval '(cd $xdir && ${LN_S} $baseobj $oldobj)' || exit $?
	  fi
	done

	# Use standard objects if they are pic
	test -z "$pic_flag" && libobjs=`$echo "X$libobjs" | $SP2NL | $Xsed -e "$lo2o" | $NL2SP`

	# Prepare the list of exported symbols
	if test -z "$export_symbols"; then
	  if test "$always_export_symbols" = yes || test -n "$export_symbols_regex"; then
	    $show "generating symbol list for \`$libname.la'"
	    export_symbols="$output_objdir/$libname.exp"
	    $run $rm $export_symbols
	    eval cmds=\"$export_symbols_cmds\"
	    save_ifs="$IFS"; IFS='~'
	    for cmd in $cmds; do
	      IFS="$save_ifs"
	      $show "$cmd"
	      $run eval "$cmd" || exit $?
	    done
	    IFS="$save_ifs"
	    if test -n "$export_symbols_regex"; then
	      $show "egrep -e \"$export_symbols_regex\" \"$export_symbols\" > \"${export_symbols}T\""
	      $run eval 'egrep -e "$export_symbols_regex" "$export_symbols" > "${export_symbols}T"'
	      $show "$mv \"${export_symbols}T\" \"$export_symbols\""
	      $run eval '$mv "${export_symbols}T" "$export_symbols"'
	    fi
	  fi
	fi

	if test -n "$export_symbols" && test -n "$include_expsyms"; then
	  $run eval '$echo "X$include_expsyms" | $SP2NL >> "$export_symbols"'
	fi

	if test -n "$convenience"; then
	  if test -n "$whole_archive_flag_spec"; then
	    eval libobjs=\"\$libobjs $whole_archive_flag_spec\"
	  else
	    gentop="$output_objdir/${outputname}x"
	    $show "${rm}r $gentop"
	    $run ${rm}r "$gentop"
	    $show "mkdir $gentop"
	    $run mkdir "$gentop"
	    status=$?
	    if test $status -ne 0 && test ! -d "$gentop"; then
	      exit $status
	    fi
	    generated="$generated $gentop"

	    for xlib in $convenience; do
	      # Extract the objects.
	      case $xlib in
	      [\\/]* | [A-Za-z]:[\\/]*) xabs="$xlib" ;;
	      *) xabs=`pwd`"/$xlib" ;;
	      esac
	      xlib=`$echo "X$xlib" | $Xsed -e 's%^.*/%%'`
	      xdir="$gentop/$xlib"

	      $show "${rm}r $xdir"
	      $run ${rm}r "$xdir"
	      $show "mkdir $xdir"
	      $run mkdir "$xdir"
	      status=$?
	      if test $status -ne 0 && test ! -d "$xdir"; then
		exit $status
	      fi
	      $show "(cd $xdir && $AR x $xabs)"
	      $run eval "(cd \$xdir && $AR x \$xabs)" || exit $?

	      libobjs="$libobjs "`find $xdir -name \*.o -print -o -name \*.lo -print | $NL2SP`
	    done
	  fi
	fi

	if test "$thread_safe" = yes && test -n "$thread_safe_flag_spec"; then
	  eval flag=\"$thread_safe_flag_spec\"
	  linker_flags="$linker_flags $flag"
	fi

	# Make a backup of the uninstalled library when relinking
	if test "$mode" = relink; then
	  $run eval '(cd $output_objdir && $rm ${realname}U && $mv $realname ${realname}U)' || exit $?
	fi

	# Do each of the archive commands.
	if test -n "$export_symbols" && test -n "$archive_expsym_cmds"; then
	  eval cmds=\"$archive_expsym_cmds\"
	else
	  save_deplibs="$deplibs"
	  for conv in $convenience; do
	    tmp_deplibs=
	    for test_deplib in $deplibs; do
	      if test "$test_deplib" != "$conv"; then
		tmp_deplibs="$tmp_deplibs $test_deplib"
	      fi
	    done
	    deplibs="$tmp_deplibs"
	  done
	  eval cmds=\"$archive_cmds\"
	  deplibs="$save_deplibs"
	fi
	save_ifs="$IFS"; IFS='~'
	for cmd in $cmds; do
	  IFS="$save_ifs"
	  $show "$cmd"
	  $run eval "$cmd" || exit $?
	done
	IFS="$save_ifs"

	# Restore the uninstalled library and exit
	if test "$mode" = relink; then
	  $run eval '(cd $output_objdir && $rm ${realname}T && $mv $realname ${realname}T && $mv "$realname"U $realname)' || exit $?
	  exit 0
	fi

	# Create links to the real library.
	for linkname in $linknames; do
	  if test "$realname" != "$linkname"; then
	    $show "(cd $output_objdir && $rm $linkname && $LN_S $realname $linkname)"
	    $run eval '(cd $output_objdir && $rm $linkname && $LN_S $realname $linkname)' || exit $?
	  fi
	done

	# If -module or -export-dynamic was specified, set the dlname.
	if test "$module" = yes || test "$export_dynamic" = yes; then
	  # On all known operating systems, these are identical.
	  dlname="$soname"
	fi
      fi
      ;;

    obj)
      if test -n "$deplibs"; then
	$echo "$modename: warning: \`-l' and \`-L' are ignored for objects" 1>&2
      fi

      if test -n "$dlfiles$dlprefiles" || test "$dlself" != no; then
	$echo "$modename: warning: \`-dlopen' is ignored for objects" 1>&2
      fi

      if test -n "$rpath"; then
	$echo "$modename: warning: \`-rpath' is ignored for objects" 1>&2
      fi

      if test -n "$xrpath"; then
	$echo "$modename: warning: \`-R' is ignored for objects" 1>&2
      fi

      if test -n "$vinfo"; then
	$echo "$modename: warning: \`-version-info' is ignored for objects" 1>&2
      fi

      if test -n "$release"; then
	$echo "$modename: warning: \`-release' is ignored for objects" 1>&2
      fi

      case $output in
      *.lo)
	if test -n "$objs$old_deplibs"; then
	  $echo "$modename: cannot build library object \`$output' from non-libtool objects" 1>&2
	  exit 1
	fi
	libobj="$output"
	obj=`$echo "X$output" | $Xsed -e "$lo2o"`
	;;
      *)
	libobj=
	obj="$output"
	;;
      esac

      # Delete the old objects.
      $run $rm $obj $libobj

      # Objects from convenience libraries.  This assumes
      # single-version convenience libraries.  Whenever we create
      # different ones for PIC/non-PIC, this we'll have to duplicate
      # the extraction.
      reload_conv_objs=
      gentop=
      # reload_cmds runs $LD directly, so let us get rid of
      # -Wl from whole_archive_flag_spec
      wl=

      if test -n "$convenience"; then
	if test -n "$whole_archive_flag_spec"; then
	  eval reload_conv_objs=\"\$reload_objs $whole_archive_flag_spec\"
	else
	  gentop="$output_objdir/${obj}x"
	  $show "${rm}r $gentop"
	  $run ${rm}r "$gentop"
	  $show "mkdir $gentop"
	  $run mkdir "$gentop"
	  status=$?
	  if test $status -ne 0 && test ! -d "$gentop"; then
	    exit $status
	  fi
	  generated="$generated $gentop"

	  for xlib in $convenience; do
	    # Extract the objects.
	    case $xlib in
	    [\\/]* | [A-Za-z]:[\\/]*) xabs="$xlib" ;;
	    *) xabs=`pwd`"/$xlib" ;;
	    esac
	    xlib=`$echo "X$xlib" | $Xsed -e 's%^.*/%%'`
	    xdir="$gentop/$xlib"

	    $show "${rm}r $xdir"
	    $run ${rm}r "$xdir"
	    $show "mkdir $xdir"
	    $run mkdir "$xdir"
	    status=$?
	    if test $status -ne 0 && test ! -d "$xdir"; then
	      exit $status
	    fi
	    $show "(cd $xdir && $AR x $xabs)"
	    $run eval "(cd \$xdir && $AR x \$xabs)" || exit $?

	    reload_conv_objs="$reload_objs "`find $xdir -name \*.o -print -o -name \*.lo -print | $NL2SP`
	  done
	fi
      fi

      # Create the old-style object.
      reload_objs="$objs$old_deplibs "`$echo "X$libobjs" | $SP2NL | $Xsed -e '/\.'${libext}$'/d' -e '/\.lib$/d' -e "$lo2o" | $NL2SP`" $reload_conv_objs" ### testsuite: skip nested quoting test

      output="$obj"
      eval cmds=\"$reload_cmds\"
      save_ifs="$IFS"; IFS='~'
      for cmd in $cmds; do
	IFS="$save_ifs"
	$show "$cmd"
	$run eval "$cmd" || exit $?
      done
      IFS="$save_ifs"

      # Exit if we aren't doing a library object file.
      if test -z "$libobj"; then
	if test -n "$gentop"; then
	  $show "${rm}r $gentop"
	  $run ${rm}r $gentop
	fi

	exit 0
      fi

      if test "$build_libtool_libs" != yes; then
	if test -n "$gentop"; then
	  $show "${rm}r $gentop"
	  $run ${rm}r $gentop
	fi

	# Create an invalid libtool object if no PIC, so that we don't
	# accidentally link it into a program.
	$show "echo timestamp > $libobj"
	$run eval "echo timestamp > $libobj" || exit $?
	exit 0
      fi

      if test -n "$pic_flag" || test "$pic_mode" != default; then
	# Only do commands if we really have different PIC objects.
	reload_objs="$libobjs $reload_conv_objs"
	output="$libobj"
	eval cmds=\"$reload_cmds\"
	save_ifs="$IFS"; IFS='~'
	for cmd in $cmds; do
	  IFS="$save_ifs"
	  $show "$cmd"
	  $run eval "$cmd" || exit $?
	done
	IFS="$save_ifs"
      else
	# Just create a symlink.
	$show $rm $libobj
	$run $rm $libobj
	xdir=`$echo "X$libobj" | $Xsed -e 's%/[^/]*$%%'`
	if test "X$xdir" = "X$libobj"; then
	  xdir="."
	else
	  xdir="$xdir"
	fi
	baseobj=`$echo "X$libobj" | $Xsed -e 's%^.*/%%'`
	oldobj=`$echo "X$baseobj" | $Xsed -e "$lo2o"`
	$show "(cd $xdir && $LN_S $oldobj $baseobj)"
	$run eval '(cd $xdir && $LN_S $oldobj $baseobj)' || exit $?
      fi

      if test -n "$gentop"; then
	$show "${rm}r $gentop"
	$run ${rm}r $gentop
      fi

      exit 0
      ;;

    prog)
      case $host in
	*cygwin*) output=`echo $output | ${SED} -e 's,.exe$,,;s,$,.exe,'` ;;
      esac
      if test -n "$vinfo"; then
	$echo "$modename: warning: \`-version-info' is ignored for programs" 1>&2
      fi

      if test -n "$release"; then
	$echo "$modename: warning: \`-release' is ignored for programs" 1>&2
      fi

      if test "$preload" = yes; then
	if test "$dlopen_support" = unknown && test "$dlopen_self" = unknown &&
	   test "$dlopen_self_static" = unknown; then
	  $echo "$modename: warning: \`AC_LIBTOOL_DLOPEN' not used. Assuming no dlopen support."
	fi
      fi

      case $host in
      *-*-rhapsody* | *-*-darwin1.[012])
	# On Rhapsody replace the C library is the System framework
	compile_deplibs=`$echo "X $compile_deplibs" | $Xsed -e 's/ -lc / -framework System /'`
	finalize_deplibs=`$echo "X $finalize_deplibs" | $Xsed -e 's/ -lc / -framework System /'`
	case $host in
	*darwin*)
	  # Don't allow lazy linking, it breaks C++ global constructors
	  compile_command="$compile_command ${wl}-bind_at_load"
	  finalize_command="$finalize_command ${wl}-bind_at_load"
	  ;;
	esac
	;;
      esac

      compile_command="$compile_command $compile_deplibs"
      finalize_command="$finalize_command $finalize_deplibs"

      if test -n "$rpath$xrpath"; then
	# If the user specified any rpath flags, then add them.
	for libdir in $rpath $xrpath; do
	  # This is the magic to use -rpath.
	  case "$finalize_rpath " in
	  *" $libdir "*) ;;
	  *) finalize_rpath="$finalize_rpath $libdir" ;;
	  esac
	done
      fi

      # Now hardcode the library paths
      rpath=
      hardcode_libdirs=
      for libdir in $compile_rpath $finalize_rpath; do
	if test -n "$hardcode_libdir_flag_spec"; then
	  if test -n "$hardcode_libdir_separator"; then
	    if test -z "$hardcode_libdirs"; then
	      hardcode_libdirs="$libdir"
	    else
	      # Just accumulate the unique libdirs.
	      case $hardcode_libdir_separator$hardcode_libdirs$hardcode_libdir_separator in
	      *"$hardcode_libdir_separator$libdir$hardcode_libdir_separator"*)
		;;
	      *)
		hardcode_libdirs="$hardcode_libdirs$hardcode_libdir_separator$libdir"
		;;
	      esac
	    fi
	  else
	    eval flag=\"$hardcode_libdir_flag_spec\"
	    rpath="$rpath $flag"
	  fi
	elif test -n "$runpath_var"; then
	  case "$perm_rpath " in
	  *" $libdir "*) ;;
	  *) perm_rpath="$perm_rpath $libdir" ;;
	  esac
	fi
	case $host in
	*-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
	  case :$dllsearchpath: in
	  *":$libdir:"*) ;;
	  *) dllsearchpath="$dllsearchpath:$libdir";;
	  esac
	  ;;
	esac
      done
      # Substitute the hardcoded libdirs into the rpath.
      if test -n "$hardcode_libdir_separator" &&
	 test -n "$hardcode_libdirs"; then
	libdir="$hardcode_libdirs"
	eval rpath=\" $hardcode_libdir_flag_spec\"
      fi
      compile_rpath="$rpath"

      rpath=
      hardcode_libdirs=
      for libdir in $finalize_rpath; do
	if test -n "$hardcode_libdir_flag_spec"; then
	  if test -n "$hardcode_libdir_separator"; then
	    if test -z "$hardcode_libdirs"; then
	      hardcode_libdirs="$libdir"
	    else
	      # Just accumulate the unique libdirs.
	      case $hardcode_libdir_separator$hardcode_libdirs$hardcode_libdir_separator in
	      *"$hardcode_libdir_separator$libdir$hardcode_libdir_separator"*)
		;;
	      *)
		hardcode_libdirs="$hardcode_libdirs$hardcode_libdir_separator$libdir"
		;;
	      esac
	    fi
	  else
	    eval flag=\"$hardcode_libdir_flag_spec\"
	    rpath="$rpath $flag"
	  fi
	elif test -n "$runpath_var"; then
	  case "$finalize_perm_rpath " in
	  *" $libdir "*) ;;
	  *) finalize_perm_rpath="$finalize_perm_rpath $libdir" ;;
	  esac
	fi
      done
      # Substitute the hardcoded libdirs into the rpath.
      if test -n "$hardcode_libdir_separator" &&
	 test -n "$hardcode_libdirs"; then
	libdir="$hardcode_libdirs"
	eval rpath=\" $hardcode_libdir_flag_spec\"
      fi
      finalize_rpath="$rpath"

      if test -n "$libobjs" && test "$build_old_libs" = yes; then
	# Transform all the library objects into standard objects.
	compile_command=`$echo "X$compile_command" | $SP2NL | $Xsed -e "$lo2o" | $NL2SP`
	finalize_command=`$echo "X$finalize_command" | $SP2NL | $Xsed -e "$lo2o" | $NL2SP`
      fi

      dlsyms=
      if test -n "$dlfiles$dlprefiles" || test "$dlself" != no; then
	if test -n "$NM" && test -n "$global_symbol_pipe"; then
	  dlsyms="${outputname}S.c"
	else
	  $echo "$modename: not configured to extract global symbols from dlpreopened files" 1>&2
	fi
      fi

      if test -n "$dlsyms"; then
	case $dlsyms in
	"") ;;
	*.c)
	  # Discover the nlist of each of the dlfiles.
	  nlist="$output_objdir/${outputname}.nm"

	  $show "$rm $nlist ${nlist}S ${nlist}T"
	  $run $rm "$nlist" "${nlist}S" "${nlist}T"

	  # Parse the name list into a source file.
	  $show "creating $output_objdir/$dlsyms"

	  test -z "$run" && $echo > "$output_objdir/$dlsyms" "\
/* $dlsyms - symbol resolution table for \`$outputname' dlsym emulation. */
/* Generated by $PROGRAM - GNU $PACKAGE $VERSION$TIMESTAMP */

#ifdef __cplusplus
extern \"C\" {
#endif

/* Prevent the only kind of declaration conflicts we can make. */
#define lt_preloaded_symbols some_other_symbol

/* External symbol declarations for the compiler. */\
"

	  if test "$dlself" = yes; then
	    $show "generating symbol list for \`$output'"

	    test -z "$run" && $echo ': @PROGRAM@ ' > "$nlist"

	    # Add our own program objects to the symbol list.
	    progfiles=`$echo "X$objs$old_deplibs" | $SP2NL | $Xsed -e "$lo2o" | $NL2SP`
	    for arg in $progfiles; do
	      $show "extracting global C symbols from \`$arg'"
	      $run eval "$NM $arg | $global_symbol_pipe >> '$nlist'"
	    done

	    if test -n "$exclude_expsyms"; then
	      $run eval 'egrep -v " ($exclude_expsyms)$" "$nlist" > "$nlist"T'
	      $run eval '$mv "$nlist"T "$nlist"'
	    fi

	    if test -n "$export_symbols_regex"; then
	      $run eval 'egrep -e "$export_symbols_regex" "$nlist" > "$nlist"T'
	      $run eval '$mv "$nlist"T "$nlist"'
	    fi

	    # Prepare the list of exported symbols
	    if test -z "$export_symbols"; then
	      export_symbols="$output_objdir/$output.exp"
	      $run $rm $export_symbols
	      $run eval "${SED} -n -e '/^: @PROGRAM@$/d' -e 's/^.* \(.*\)$/\1/p' "'< "$nlist" > "$export_symbols"'
	    else
	      $run eval "${SED} -e 's/\([][.*^$]\)/\\\1/g' -e 's/^/ /' -e 's/$/$/'"' < "$export_symbols" > "$output_objdir/$output.exp"'
	      $run eval 'grep -f "$output_objdir/$output.exp" < "$nlist" > "$nlist"T'
	      $run eval 'mv "$nlist"T "$nlist"'
	    fi
	  fi

	  for arg in $dlprefiles; do
	    $show "extracting global C symbols from \`$arg'"
	    name=`echo "$arg" | ${SED} -e 's%^.*/%%'`
	    $run eval 'echo ": $name " >> "$nlist"'
	    $run eval "$NM $arg | $global_symbol_pipe >> '$nlist'"
	  done

	  if test -z "$run"; then
	    # Make sure we have at least an empty file.
	    test -f "$nlist" || : > "$nlist"

	    if test -n "$exclude_expsyms"; then
	      egrep -v " ($exclude_expsyms)$" "$nlist" > "$nlist"T
	      $mv "$nlist"T "$nlist"
	    fi

	    # Try sorting and uniquifying the output.
	    if grep -v "^: " < "$nlist" |
		if sort -k 3 </dev/null >/dev/null 2>&1; then
		  sort -k 3
		else
		  sort +2
		fi |
		uniq > "$nlist"S; then
	      :
	    else
	      grep -v "^: " < "$nlist" > "$nlist"S
	    fi

	    if test -f "$nlist"S; then
	      eval "$global_symbol_to_cdecl"' < "$nlist"S >> "$output_objdir/$dlsyms"'
	    else
	      echo '/* NONE */' >> "$output_objdir/$dlsyms"
	    fi

	    $echo >> "$output_objdir/$dlsyms" "\

#undef lt_preloaded_symbols

#if defined (__STDC__) && __STDC__
# define lt_ptr void *
#else
# define lt_ptr char *
# define const
#endif

/* The mapping between symbol names and symbols. */
const struct {
  const char *name;
  lt_ptr address;
}
lt_preloaded_symbols[] =
{\
"

	    eval "$global_symbol_to_c_name_address" < "$nlist" >> "$output_objdir/$dlsyms"

	    $echo >> "$output_objdir/$dlsyms" "\
  {0, (lt_ptr) 0}
};

/* This works around a problem in FreeBSD linker */
#ifdef FREEBSD_WORKAROUND
static const void *lt_preloaded_setup() {
  return lt_preloaded_symbols;
}
#endif

#ifdef __cplusplus
}
#endif\
"
	  fi

	  pic_flag_for_symtable=
	  case $host in
	  # compiling the symbol table file with pic_flag works around
	  # a FreeBSD bug that causes programs to crash when -lm is
	  # linked before any other PIC object.  But we must not use
	  # pic_flag when linking with -static.  The problem exists in
	  # FreeBSD 2.2.6 and is fixed in FreeBSD 3.1.
	  *-*-freebsd2*|*-*-freebsd3.0*|*-*-freebsdelf3.0*)
	    case "$compile_command " in
	    *" -static "*) ;;
	    *) pic_flag_for_symtable=" $pic_flag -DPIC -DFREEBSD_WORKAROUND";;
	    esac;;
	  *-*-hpux*)
	    case "$compile_command " in
	    *" -static "*) ;;
	    *) pic_flag_for_symtable=" $pic_flag -DPIC";;
	    esac
	  esac

	  # Now compile the dynamic symbol file.
	  $show "(cd $output_objdir && $CC -c$no_builtin_flag$pic_flag_for_symtable \"$dlsyms\")"
	  $run eval '(cd $output_objdir && $CC -c$no_builtin_flag$pic_flag_for_symtable "$dlsyms")' || exit $?

	  # Clean up the generated files.
	  $show "$rm $output_objdir/$dlsyms $nlist ${nlist}S ${nlist}T"
	  $run $rm "$output_objdir/$dlsyms" "$nlist" "${nlist}S" "${nlist}T"

	  # Transform the symbol file into the correct name.
	  compile_command=`$echo "X$compile_command" | $Xsed -e "s%@SYMFILE@%$output_objdir/${outputname}S.${objext}%"`
	  finalize_command=`$echo "X$finalize_command" | $Xsed -e "s%@SYMFILE@%$output_objdir/${outputname}S.${objext}%"`
	  ;;
	*)
	  $echo "$modename: unknown suffix for \`$dlsyms'" 1>&2
	  exit 1
	  ;;
	esac
      else
	# We keep going just in case the user didn't refer to
	# lt_preloaded_symbols.  The linker will fail if global_symbol_pipe
	# really was required.

	# Nullify the symbol file.
	compile_command=`$echo "X$compile_command" | $Xsed -e "s% @SYMFILE@%%"`
	finalize_command=`$echo "X$finalize_command" | $Xsed -e "s% @SYMFILE@%%"`
      fi

      if test $need_relink = no || test "$build_libtool_libs" != yes; then
	# Replace the output file specification.
	compile_command=`$echo "X$compile_command" | $Xsed -e 's%@OUTPUT@%'"$output"'%g'`
	link_command="$compile_command$compile_rpath"

	# We have no uninstalled library dependencies, so finalize right now.
	$show "$link_command"
	$run eval "$link_command"
	status=$?

	# Delete the generated files.
	if test -n "$dlsyms"; then
	  $show "$rm $output_objdir/${outputname}S.${objext}"
	  $run $rm "$output_objdir/${outputname}S.${objext}"
	fi

	exit $status
      fi

      if test -n "$shlibpath_var"; then
	# We should set the shlibpath_var
	rpath=
	for dir in $temp_rpath; do
	  case $dir in
	  [\\/]* | [A-Za-z]:[\\/]*)
	    # Absolute path.
	    rpath="$rpath$dir:"
	    ;;
	  *)
	    # Relative path: add a thisdir entry.
	    rpath="$rpath\$thisdir/$dir:"
	    ;;
	  esac
	done
	temp_rpath="$rpath"
      fi

      if test -n "$compile_shlibpath$finalize_shlibpath"; then
	compile_command="$shlibpath_var=\"$compile_shlibpath$finalize_shlibpath\$$shlibpath_var\" $compile_command"
      fi
      if test -n "$finalize_shlibpath"; then
	finalize_command="$shlibpath_var=\"$finalize_shlibpath\$$shlibpath_var\" $finalize_command"
      fi

      compile_var=
      finalize_var=
      if test -n "$runpath_var"; then
	if test -n "$perm_rpath"; then
	  # We should set the runpath_var.
	  rpath=
	  for dir in $perm_rpath; do
	    rpath="$rpath$dir:"
	  done
	  compile_var="$runpath_var=\"$rpath\$$runpath_var\" "
	fi
	if test -n "$finalize_perm_rpath"; then
	  # We should set the runpath_var.
	  rpath=
	  for dir in $finalize_perm_rpath; do
	    rpath="$rpath$dir:"
	  done
	  finalize_var="$runpath_var=\"$rpath\$$runpath_var\" "
	fi
      fi

      if test "$no_install" = yes; then
	# We don't need to create a wrapper script.
	link_command="$compile_var$compile_command$compile_rpath"
	# Replace the output file specification.
	link_command=`$echo "X$link_command" | $Xsed -e 's%@OUTPUT@%'"$output"'%g'`
	# Delete the old output file.
	$run $rm $output
	# Link the executable and exit
	$show "$link_command"
	$run eval "$link_command" || exit $?
	exit 0
      fi

      if test "$hardcode_action" = relink; then
	# Fast installation is not supported
	link_command="$compile_var$compile_command$compile_rpath"
	relink_command="$finalize_var$finalize_command$finalize_rpath"

	$echo "$modename: warning: this platform does not like uninstalled shared libraries" 1>&2
	$echo "$modename: \`$output' will be relinked during installation" 1>&2
      else
	if test "$fast_install" != no; then
	  link_command="$finalize_var$compile_command$finalize_rpath"
	  if test "$fast_install" = yes; then
	    relink_command=`$echo "X$compile_var$compile_command$compile_rpath" | $Xsed -e 's%@OUTPUT@%\$progdir/\$file%g'`
	  else
	    # fast_install is set to needless
	    relink_command=
	  fi
	else
	  link_command="$compile_var$compile_command$compile_rpath"
	  relink_command="$finalize_var$finalize_command$finalize_rpath"
	fi
      fi

      # Replace the output file specification.
      link_command=`$echo "X$link_command" | $Xsed -e 's%@OUTPUT@%'"$output_objdir/$outputname"'%g'`

      # Delete the old output files.
      $run $rm $output $output_objdir/$outputname $output_objdir/lt-$outputname

      $show "$link_command"
      $run eval "$link_command" || exit $?

      # Now create the wrapper script.
      $show "creating $output"

      # Quote the relink command for shipping.
      if test -n "$relink_command"; then
	# Preserve any variables that may affect compiler behavior
	for var in $variables_saved_for_relink; do
	  if eval test -z \"\${$var+set}\"; then
	    relink_command="{ test -z \"\${$var+set}\" || unset $var || { $var=; export $var; }; }; $relink_command"
	  elif eval var_value=\$$var; test -z "$var_value"; then
	    relink_command="$var=; export $var; $relink_command"
	  else
	    var_value=`$echo "X$var_value" | $Xsed -e "$sed_quote_subst"`
	    relink_command="$var=\"$var_value\"; export $var; $relink_command"
	  fi
	done
	relink_command="(cd `pwd`; $relink_command)"
	relink_command=`$echo "X$relink_command" | $Xsed -e "$sed_quote_subst"`
      fi

      # Quote $echo for shipping.
      if test "X$echo" = "X$SHELL $0 --fallback-echo"; then
	case $0 in
	[\\/]* | [A-Za-z]:[\\/]*) qecho="$SHELL $0 --fallback-echo";;
	*) qecho="$SHELL `pwd`/$0 --fallback-echo";;
	esac
	qecho=`$echo "X$qecho" | $Xsed -e "$sed_quote_subst"`
      else
	qecho=`$echo "X$echo" | $Xsed -e "$sed_quote_subst"`
      fi

      # Only actually do things if our run command is non-null.
      if test -z "$run"; then
	# win32 will think the script is a binary if it has
	# a .exe suffix, so we strip it off here.
	case $output in
	  *.exe) output=`echo $output|${SED} 's,.exe$,,'` ;;
	esac
	# test for cygwin because mv fails w/o .exe extensions
	case $host in
	  *cygwin*) exeext=.exe ;;
	  *) exeext= ;;
	esac
	$rm $output
	trap "$rm $output; exit 1" 1 2 15

	$echo > $output "\
#! $SHELL

# $output - temporary wrapper script for $objdir/$outputname
# Generated by $PROGRAM - GNU $PACKAGE $VERSION$TIMESTAMP
#
# The $output program cannot be directly executed until all the libtool
# libraries that it depends on are installed.
#
# This wrapper script should never be moved out of the build directory.
# If it is, it will not operate correctly.

# Sed substitution that helps us do robust quoting.  It backslashifies
# metacharacters that are still active within double-quoted strings.
Xsed="${SED}"' -e 1s/^X//'
sed_quote_subst='$sed_quote_subst'

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
if test \"\${CDPATH+set}\" = set; then CDPATH=:; export CDPATH; fi

relink_command=\"$relink_command\"

# This environment variable determines our operation mode.
if test \"\$libtool_install_magic\" = \"$magic\"; then
  # install mode needs the following variable:
  notinst_deplibs='$notinst_deplibs'
else
  # When we are sourced in execute mode, \$file and \$echo are already set.
  if test \"\$libtool_execute_magic\" != \"$magic\"; then
    echo=\"$qecho\"
    file=\"\$0\"
    # Make sure echo works.
    if test \"X\$1\" = X--no-reexec; then
      # Discard the --no-reexec flag, and continue.
      shift
    elif test \"X\`(\$echo '\t') 2>/dev/null\`\" = 'X\t'; then
      # Yippee, \$echo works!
      :
    else
      # Restart under the correct shell, and then maybe \$echo will work.
      exec $SHELL \"\$0\" --no-reexec \${1+\"\$@\"}
    fi
  fi\
"
	$echo >> $output "\

  # Find the directory that this script lives in.
  thisdir=\`\$echo \"X\$file\" | \$Xsed -e 's%/[^/]*$%%'\`
  test \"x\$thisdir\" = \"x\$file\" && thisdir=.

  # Follow symbolic links until we get to the real thisdir.
  file=\`ls -ld \"\$file\" | ${SED} -n 's/.*-> //p'\`
  while test -n \"\$file\"; do
    destdir=\`\$echo \"X\$file\" | \$Xsed -e 's%/[^/]*\$%%'\`

    # If there was a directory component, then change thisdir.
    if test \"x\$destdir\" != \"x\$file\"; then
      case \"\$destdir\" in
      [\\\\/]* | [A-Za-z]:[\\\\/]*) thisdir=\"\$destdir\" ;;
      *) thisdir=\"\$thisdir/\$destdir\" ;;
      esac
    fi

    file=\`\$echo \"X\$file\" | \$Xsed -e 's%^.*/%%'\`
    file=\`ls -ld \"\$thisdir/\$file\" | ${SED} -n 's/.*-> //p'\`
  done

  # Try to get the absolute directory name.
  absdir=\`cd \"\$thisdir\" && pwd\`
  test -n \"\$absdir\" && thisdir=\"\$absdir\"
"

	if test "$fast_install" = yes; then
	  echo >> $output "\
  program=lt-'$outputname'$exeext
  progdir=\"\$thisdir/$objdir\"

  if test ! -f \"\$progdir/\$program\" || \\
     { file=\`ls -1dt \"\$progdir/\$program\" \"\$progdir/../\$program\" 2>/dev/null | ${SED} 1q\`; \\
       test \"X\$file\" != \"X\$progdir/\$program\"; }; then

    file=\"\$\$-\$program\"

    if test ! -d \"\$progdir\"; then
      $mkdir \"\$progdir\"
    else
      $rm \"\$progdir/\$file\"
    fi"

	  echo >> $output "\

    # relink executable if necessary
    if test -n \"\$relink_command\"; then
      if relink_command_output=\`eval \$relink_command 2>&1\`; then :
      else
	$echo \"\$relink_command_output\" >&2
	$rm \"\$progdir/\$file\"
	exit 1
      fi
    fi

    $mv \"\$progdir/\$file\" \"\$progdir/\$program\" 2>/dev/null ||
    { $rm \"\$progdir/\$program\";
      $mv \"\$progdir/\$file\" \"\$progdir/\$program\"; }
    $rm \"\$progdir/\$file\"
  fi"
	else
	  echo >> $output "\
  program='$outputname'
  progdir=\"\$thisdir/$objdir\"
"
	fi

	echo >> $output "\

  if test -f \"\$progdir/\$program\"; then"

	# Export our shlibpath_var if we have one.
	if test "$shlibpath_overrides_runpath" = yes && test -n "$shlibpath_var" && test -n "$temp_rpath"; then
	  $echo >> $output "\
    # Add our own library path to $shlibpath_var
    $shlibpath_var=\"$temp_rpath\$$shlibpath_var\"

    # Some systems cannot cope with colon-terminated $shlibpath_var
    # The second colon is a workaround for a bug in BeOS R4 ${SED}
    $shlibpath_var=\`\$echo \"X\$$shlibpath_var\" | \$Xsed -e 's/::*\$//'\`

    export $shlibpath_var
"
	fi

	# fixup the dll searchpath if we need to.
	if test -n "$dllsearchpath"; then
	  $echo >> $output "\
    # Add the dll search path components to the executable PATH
    PATH=$dllsearchpath:\$PATH
"
	fi

	$echo >> $output "\
    if test \"\$libtool_execute_magic\" != \"$magic\"; then
      # Run the actual program with our arguments.
"
	case $host in
	# win32 systems need to use the prog path for dll
	# lookup to work
	*-*-cygwin* | *-*-pw32*)
	  $echo >> $output "\
      exec \$progdir/\$program \${1+\"\$@\"}
"
	  ;;

	# Backslashes separate directories on plain windows
	*-*-mingw | *-*-os2*)
	  $echo >> $output "\
      exec \$progdir\\\\\$program \${1+\"\$@\"}
"
	  ;;

	*)
	  $echo >> $output "\
      # Export the path to the program.
      PATH=\"\$progdir:\$PATH\"
      export PATH

      exec \$program \${1+\"\$@\"}
"
	  ;;
	esac
	$echo >> $output "\
      \$echo \"\$0: cannot exec \$program \${1+\"\$@\"}\"
      exit 1
    fi
  else
    # The program doesn't exist.
    \$echo \"\$0: error: \$progdir/\$program does not exist\" 1>&2
    \$echo \"This script is just a wrapper for \$program.\" 1>&2
    echo \"See the $PACKAGE documentation for more information.\" 1>&2
    exit 1
  fi
fi\
"
	chmod +x $output
      fi
      exit 0
      ;;
    esac

    # See if we need to build an old-fashioned archive.
    for oldlib in $oldlibs; do

      if test "$build_libtool_libs" = convenience; then
	oldobjs="$libobjs_save"
	addlibs="$convenience"
	build_libtool_libs=no
      else
	if test "$build_libtool_libs" = module; then
	  oldobjs="$libobjs_save"
	  build_libtool_libs=no
	else
	  oldobjs="$objs$old_deplibs "`$echo "X$libobjs_save" | $SP2NL | $Xsed -e '/\.'${libext}'$/d' -e '/\.lib$/d' -e "$lo2o" | $NL2SP`
	fi
	addlibs="$old_convenience"
      fi

      if test -n "$addlibs"; then
	gentop="$output_objdir/${outputname}x"
	$show "${rm}r $gentop"
	$run ${rm}r "$gentop"
	$show "mkdir $gentop"
	$run mkdir "$gentop"
	status=$?
	if test $status -ne 0 && test ! -d "$gentop"; then
	  exit $status
	fi
	generated="$generated $gentop"

	# Add in members from convenience archives.
	for xlib in $addlibs; do
	  # Extract the objects.
	  case $xlib in
	  [\\/]* | [A-Za-z]:[\\/]*) xabs="$xlib" ;;
	  *) xabs=`pwd`"/$xlib" ;;
	  esac
	  xlib=`$echo "X$xlib" | $Xsed -e 's%^.*/%%'`
	  xdir="$gentop/$xlib"

	  $show "${rm}r $xdir"
	  $run ${rm}r "$xdir"
	  $show "mkdir $xdir"
	  $run mkdir "$xdir"
	  status=$?
	  if test $status -ne 0 && test ! -d "$xdir"; then
	    exit $status
	  fi
	  $show "(cd $xdir && $AR x $xabs)"
	  $run eval "(cd \$xdir && $AR x \$xabs)" || exit $?

	  oldobjs="$oldobjs "`find $xdir -name \*.${objext} -print -o -name \*.lo -print | $NL2SP`
	done
      fi

      # Do each command in the archive commands.
      if test -n "$old_archive_from_new_cmds" && test "$build_libtool_libs" = yes; then
	eval cmds=\"$old_archive_from_new_cmds\"
      else
	# Ensure that we have .o objects in place in case we decided
	# not to build a shared library, and have fallen back to building
	# static libs even though --disable-static was passed!
	for oldobj in $oldobjs; do
	  if test ! -f $oldobj; then
	    xdir=`$echo "X$oldobj" | $Xsed -e 's%/[^/]*$%%'`
	    if test "X$xdir" = "X$oldobj"; then
	      xdir="."
	    else
	      xdir="$xdir"
	    fi
	    baseobj=`$echo "X$oldobj" | $Xsed -e 's%^.*/%%'`
	    obj=`$echo "X$baseobj" | $Xsed -e "$o2lo"`
	    $show "(cd $xdir && ${LN_S} $obj $baseobj)"
	    $run eval '(cd $xdir && ${LN_S} $obj $baseobj)' || exit $?
	  fi
	done

	eval cmds=\"$old_archive_cmds\"
      fi
      save_ifs="$IFS"; IFS='~'
      for cmd in $cmds; do
	IFS="$save_ifs"
	$show "$cmd"
	$run eval "$cmd" || exit $?
      done
      IFS="$save_ifs"
    done

    if test -n "$generated"; then
      $show "${rm}r$generated"
      $run ${rm}r$generated
    fi

    # Now create the libtool archive.
    case $output in
    *.la)
      old_library=
      test "$build_old_libs" = yes && old_library="$libname.$libext"
      $show "creating $output"

      # Preserve any variables that may affect compiler behavior
      for var in $variables_saved_for_relink; do
	if eval test -z \"\${$var+set}\"; then
	  relink_command="{ test -z \"\${$var+set}\" || unset $var || { $var=; export $var; }; }; $relink_command"
	elif eval var_value=\$$var; test -z "$var_value"; then
	  relink_command="$var=; export $var; $relink_command"
	else
	  var_value=`$echo "X$var_value" | $Xsed -e "$sed_quote_subst"`
	  relink_command="$var=\"$var_value\"; export $var; $relink_command"
	fi
      done
      # Quote the link command for shipping.
      relink_command="(cd `pwd`; $SHELL $0 --mode=relink $libtool_args)"
      relink_command=`$echo "X$relink_command" | $Xsed -e "$sed_quote_subst"`

      # Only create the output if not a dry run.
      if test -z "$run"; then
	for installed in no yes; do
	  if test "$installed" = yes; then
	    if test -z "$install_libdir"; then
	      break
	    fi
	    output="$output_objdir/$outputname"i
	    # Replace all uninstalled libtool libraries with the installed ones
	    newdependency_libs=
	    for deplib in $dependency_libs; do
	      case $deplib in
	      *.la)
		name=`$echo "X$deplib" | $Xsed -e 's%^.*/%%'`
		eval libdir=`${SED} -n -e 's/^libdir=\(.*\)$/\1/p' $deplib`
		if test -z "$libdir"; then
		  $echo "$modename: \`$deplib' is not a valid libtool archive" 1>&2
		  exit 1
		fi
		newdependency_libs="$newdependency_libs $libdir/$name"
		;;
	      *) newdependency_libs="$newdependency_libs $deplib" ;;
	      esac
	    done
	    dependency_libs="$newdependency_libs"
	    newdlfiles=
	    for lib in $dlfiles; do
	      name=`$echo "X$lib" | $Xsed -e 's%^.*/%%'`
	      eval libdir=`${SED} -n -e 's/^libdir=\(.*\)$/\1/p' $lib`
	      if test -z "$libdir"; then
		$echo "$modename: \`$lib' is not a valid libtool archive" 1>&2
		exit 1
	      fi
	      newdlfiles="$newdlfiles $libdir/$name"
	    done
	    dlfiles="$newdlfiles"
	    newdlprefiles=
	    for lib in $dlprefiles; do
	      name=`$echo "X$lib" | $Xsed -e 's%^.*/%%'`
	      eval libdir=`${SED} -n -e 's/^libdir=\(.*\)$/\1/p' $lib`
	      if test -z "$libdir"; then
		$echo "$modename: \`$lib' is not a valid libtool archive" 1>&2
		exit 1
	      fi
	      newdlprefiles="$newdlprefiles $libdir/$name"
	    done
	    dlprefiles="$newdlprefiles"
	  fi
	  $rm $output
	  # place dlname in correct position for cygwin
	  tdlname=$dlname
	  case $host,$output,$installed,$module,$dlname in
	    *cygwin*,*lai,yes,no,*.dll) tdlname=../bin/$dlname ;;
	  esac
	  $echo > $output "\
# $outputname - a libtool library file
# Generated by $PROGRAM - GNU $PACKAGE $VERSION$TIMESTAMP
#
# Please DO NOT delete this file!
# It is necessary for linking the library.

# The name that we can dlopen(3).
dlname='$tdlname'

# Names of this library.
library_names='$library_names'

# The name of the static archive.
old_library='$old_library'

# Libraries that this one depends upon.
dependency_libs='$dependency_libs'

# Version information for $libname.
current=$current
age=$age
revision=$revision

# Is this an already installed library?
installed=$installed

# Files to dlopen/dlpreopen
dlopen='$dlfiles'
dlpreopen='$dlprefiles'

# Directory that this library needs to be installed in:
libdir='$install_libdir'"
	  if test "$installed" = no && test $need_relink = yes; then
	    $echo >> $output "\
relink_command=\"$relink_command\""
	  fi
	done
      fi

      # Do a symbolic link so that the libtool archive can be found in
      # LD_LIBRARY_PATH before the program is installed.
      $show "(cd $output_objdir && $rm $outputname && $LN_S ../$outputname $outputname)"
      $run eval '(cd $output_objdir && $rm $outputname && $LN_S ../$outputname $outputname)' || exit $?
      ;;
    esac
    exit 0
    ;;

  # libtool install mode
  install)
    modename="$modename: install"

    # There may be an optional sh(1) argument at the beginning of
    # install_prog (especially on Windows NT).
    if test "$nonopt" = "$SHELL" || test "$nonopt" = /bin/sh ||
       # Allow the use of GNU shtool's install command.
       $echo "X$nonopt" | $Xsed | grep shtool > /dev/null; then
      # Aesthetically quote it.
      arg=`$echo "X$nonopt" | $Xsed -e "$sed_quote_subst"`
      case $arg in
      *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*)
	arg="\"$arg\""
	;;
      esac
      install_prog="$arg "
      arg="$1"
      shift
    else
      install_prog=
      arg="$nonopt"
    fi

    # The real first argument should be the name of the installation program.
    # Aesthetically quote it.
    arg=`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`
    case $arg in
    *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*)
      arg="\"$arg\""
      ;;
    esac
    install_prog="$install_prog$arg"

    # We need to accept at least all the BSD install flags.
    dest=
    files=
    opts=
    prev=
    install_type=
    isdir=no
    stripme=
    for arg
    do
      if test -n "$dest"; then
	files="$files $dest"
	dest="$arg"
	continue
      fi

      case $arg in
      -d) isdir=yes ;;
      -f) prev="-f" ;;
      -g) prev="-g" ;;
      -m) prev="-m" ;;
      -o) prev="-o" ;;
      -s)
	stripme=" -s"
	continue
	;;
      -*) ;;

      *)
	# If the previous option needed an argument, then skip it.
	if test -n "$prev"; then
	  prev=
	else
	  dest="$arg"
	  continue
	fi
	;;
      esac

      # Aesthetically quote the argument.
      arg=`$echo "X$arg" | $Xsed -e "$sed_quote_subst"`
      case $arg in
      *[\[\~\#\^\&\*\(\)\{\}\|\;\<\>\?\'\ \	]*|*]*)
	arg="\"$arg\""
	;;
      esac
      install_prog="$install_prog $arg"
    done

    if test -z "$install_prog"; then
      $echo "$modename: you must specify an install program" 1>&2
      $echo "$help" 1>&2
      exit 1
    fi

    if test -n "$prev"; then
      $echo "$modename: the \`$prev' option requires an argument" 1>&2
      $echo "$help" 1>&2
      exit 1
    fi

    if test -z "$files"; then
      if test -z "$dest"; then
	$echo "$modename: no file or destination specified" 1>&2
      else
	$echo "$modename: you must specify a destination" 1>&2
      fi
      $echo "$help" 1>&2
      exit 1
    fi

    # Strip any trailing slash from the destination.
    dest=`$echo "X$dest" | $Xsed -e 's%/$%%'`

    # Check to see that the destination is a directory.
    test -d "$dest" && isdir=yes
    if test "$isdir" = yes; then
      destdir="$dest"
      destname=
    else
      destdir=`$echo "X$dest" | $Xsed -e 's%/[^/]*$%%'`
      test "X$destdir" = "X$dest" && destdir=.
      destname=`$echo "X$dest" | $Xsed -e 's%^.*/%%'`

      # Not a directory, so check to see that there is only one file specified.
      set dummy $files
      if test $# -gt 2; then
	$echo "$modename: \`$dest' is not a directory" 1>&2
	$echo "$help" 1>&2
	exit 1
      fi
    fi
    case $destdir in
    [\\/]* | [A-Za-z]:[\\/]*) ;;
    *)
      for file in $files; do
	case $file in
	*.lo) ;;
	*)
	  $echo "$modename: \`$destdir' must be an absolute directory name" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	  ;;
	esac
      done
      ;;
    esac

    # This variable tells wrapper scripts just to set variables rather
    # than running their programs.
    libtool_install_magic="$magic"

    staticlibs=
    future_libdirs=
    current_libdirs=
    for file in $files; do

      # Do each installation.
      case $file in
      *.$libext)
	# Do the static libraries later.
	staticlibs="$staticlibs $file"
	;;

      *.la)
	# Check to see that this really is a libtool archive.
	if (${SED} -e '2q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then :
	else
	  $echo "$modename: \`$file' is not a valid libtool archive" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	fi

	library_names=
	old_library=
	relink_command=
	# If there is no directory component, then add one.
	case $file in
	*/* | *\\*) . $file ;;
	*) . ./$file ;;
	esac

	# Add the libdir to current_libdirs if it is the destination.
	DESTDIR=
	if test "X$destdir" = "X$libdir"; then
	  case "$current_libdirs " in
	  *" $libdir "*) ;;
	  *) current_libdirs="$current_libdirs $libdir" ;;
	  esac
	else
	  case "$destdir" in
	    *"$libdir")
	      DESTDIR=`$echo "$destdir" | sed -e 's!'"$libdir"'$!!'`
	      if test "X$destdir" != "X$DESTDIR$libdir"; then
		DESTDIR=
	      fi
	      ;;
	  esac
	  # Note the libdir as a future libdir.
	  case "$future_libdirs " in
	  *" $libdir "*) ;;
	  *) future_libdirs="$future_libdirs $libdir" ;;
	  esac
	fi

	dir=`$echo "X$file" | $Xsed -e 's%/[^/]*$%%'`/
	test "X$dir" = "X$file/" && dir=
	dir="$dir$objdir"

	if test -n "$relink_command"; then
	  $echo "$modename: warning: relinking \`$file'" 1>&2
	  export DESTDIR
	  $show "$relink_command"
	  if $run eval "$relink_command"; then :
	  else
	    $echo "$modename: error: relink \`$file' with the above command before installing it" 1>&2
	    continue
	  fi
	fi
	unset DESTDIR

	# See the names of the shared library.
	set dummy $library_names
	if test -n "$2"; then
	  realname="$2"
	  shift
	  shift

	  srcname="$realname"
	  test -n "$relink_command" && srcname="$realname"T

	  # Install the shared library and build the symlinks.
	  $show "$install_prog $dir/$srcname $destdir/$realname"
	  $run eval "$install_prog $dir/$srcname $destdir/$realname" || exit $?
	  if test -n "$stripme" && test -n "$striplib"; then
	    $show "$striplib $destdir/$realname"
	    $run eval "$striplib $destdir/$realname" || exit $?
	  fi

	  if test $# -gt 0; then
	    # Delete the old symlinks, and create new ones.
	    for linkname
	    do
	      if test "$linkname" != "$realname"; then
		$show "(cd $destdir && $rm $linkname && $LN_S $realname $linkname)"
		$run eval "(cd $destdir && $rm $linkname && $LN_S $realname $linkname)"
	      fi
	    done
	  fi

	  # Do each command in the postinstall commands.
	  lib="$destdir/$realname"
	  eval cmds=\"$postinstall_cmds\"
	  save_ifs="$IFS"; IFS='~'
	  for cmd in $cmds; do
	    IFS="$save_ifs"
	    $show "$cmd"
	    $run eval "$cmd" || exit $?
	  done
	  IFS="$save_ifs"
	fi

	# Install the pseudo-library for information purposes.
	name=`$echo "X$file" | $Xsed -e 's%^.*/%%'`
	instname="$dir/$name"i
	$show "$install_prog $instname $destdir/$name"
	$run eval "$install_prog $instname $destdir/$name" || exit $?

	# Maybe install the static library, too.
	test -n "$old_library" && staticlibs="$staticlibs $dir/$old_library"
	;;

      *.lo)
	# Install (i.e. copy) a libtool object.

	# Figure out destination file name, if it wasn't already specified.
	if test -n "$destname"; then
	  destfile="$destdir/$destname"
	else
	  destfile=`$echo "X$file" | $Xsed -e 's%^.*/%%'`
	  destfile="$destdir/$destfile"
	fi

	# Deduce the name of the destination old-style object file.
	case $destfile in
	*.lo)
	  staticdest=`$echo "X$destfile" | $Xsed -e "$lo2o"`
	  ;;
	*.$objext)
	  staticdest="$destfile"
	  destfile=
	  ;;
	*)
	  $echo "$modename: cannot copy a libtool object to \`$destfile'" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	  ;;
	esac

	# Install the libtool object if requested.
	if test -n "$destfile"; then
	  $show "$install_prog $file $destfile"
	  $run eval "$install_prog $file $destfile" || exit $?
	fi

	# Install the old object if enabled.
	if test "$build_old_libs" = yes; then
	  # Deduce the name of the old-style object file.
	  staticobj=`$echo "X$file" | $Xsed -e "$lo2o"`

	  $show "$install_prog $staticobj $staticdest"
	  $run eval "$install_prog \$staticobj \$staticdest" || exit $?
	fi
	exit 0
	;;

      *)
	# Figure out destination file name, if it wasn't already specified.
	if test -n "$destname"; then
	  destfile="$destdir/$destname"
	else
	  destfile=`$echo "X$file" | $Xsed -e 's%^.*/%%'`
	  destfile="$destdir/$destfile"
	fi

	# Do a test to see if this is really a libtool program.
	case $host in
	*cygwin*|*mingw*)
	    wrapper=`echo $file | ${SED} -e 's,.exe$,,'`
	    ;;
	*)
	    wrapper=$file
	    ;;
	esac
	if (${SED} -e '4q' $wrapper | egrep "^# Generated by .*$PACKAGE")>/dev/null 2>&1; then
	  notinst_deplibs=
	  relink_command=

	  # If there is no directory component, then add one.
	  case $file in
	  */* | *\\*) . $wrapper ;;
	  *) . ./$wrapper ;;
	  esac

	  # Check the variables that should have been set.
	  if test -z "$notinst_deplibs"; then
	    $echo "$modename: invalid libtool wrapper script \`$wrapper'" 1>&2
	    exit 1
	  fi

	  finalize=yes
	  for lib in $notinst_deplibs; do
	    # Check to see that each library is installed.
	    libdir=
	    if test -f "$lib"; then
	      # If there is no directory component, then add one.
	      case $lib in
	      */* | *\\*) . $lib ;;
	      *) . ./$lib ;;
	      esac
	    fi
	    libfile="$libdir/"`$echo "X$lib" | $Xsed -e 's%^.*/%%g'` ### testsuite: skip nested quoting test
	    if test -n "$libdir" && test ! -f "$libfile"; then
	      $echo "$modename: warning: \`$lib' has not been installed in \`$libdir'" 1>&2
	      finalize=no
	    fi
	  done

	  relink_command=
	  # If there is no directory component, then add one.
	  case $file in
	  */* | *\\*) . $wrapper ;;
	  *) . ./$wrapper ;;
	  esac

	  outputname=
	  if test "$fast_install" = no && test -n "$relink_command"; then
	    if test "$finalize" = yes && test -z "$run"; then
	      tmpdir="/tmp"
	      test -n "$TMPDIR" && tmpdir="$TMPDIR"
              tmpdir=`mktemp -d $tmpdir/libtool-XXXXXX 2> /dev/null`
              if test $? = 0 ; then :
              else
                tmpdir="$tmpdir/libtool-$$"
              fi
	      if $mkdir -p "$tmpdir" && chmod 700 "$tmpdir"; then :
	      else
		$echo "$modename: error: cannot create temporary directory \`$tmpdir'" 1>&2
		continue
	      fi
	      file=`$echo "X$file" | $Xsed -e 's%^.*/%%'`
	      outputname="$tmpdir/$file"
	      # Replace the output file specification.
	      relink_command=`$echo "X$relink_command" | $Xsed -e 's%@OUTPUT@%'"$outputname"'%g'`

	      $show "$relink_command"
	      if $run eval "$relink_command"; then :
	      else
		$echo "$modename: error: relink \`$file' with the above command before installing it" 1>&2
		${rm}r "$tmpdir"
		continue
	      fi
	      file="$outputname"
	    else
	      $echo "$modename: warning: cannot relink \`$file'" 1>&2
	    fi
	  else
	    # Install the binary that we compiled earlier.
	    file=`$echo "X$file" | $Xsed -e "s%\([^/]*\)$%$objdir/\1%"`
	  fi
	fi

	# remove .exe since cygwin /usr/bin/install will append another
	# one anyways
	case $install_prog,$host in
	/usr/bin/install*,*cygwin*)
	  case $file:$destfile in
	  *.exe:*.exe)
	    # this is ok
	    ;;
	  *.exe:*)
	    destfile=$destfile.exe
	    ;;
	  *:*.exe)
	    destfile=`echo $destfile | ${SED} -e 's,.exe$,,'`
	    ;;
	  esac
	  ;;
	esac
	$show "$install_prog$stripme $file $destfile"
	$run eval "$install_prog\$stripme \$file \$destfile" || exit $?
	test -n "$outputname" && ${rm}r "$tmpdir"
	;;
      esac
    done

    for file in $staticlibs; do
      name=`$echo "X$file" | $Xsed -e 's%^.*/%%'`

      # Set up the ranlib parameters.
      oldlib="$destdir/$name"

      $show "$install_prog $file $oldlib"
      $run eval "$install_prog \$file \$oldlib" || exit $?

      if test -n "$stripme" && test -n "$striplib"; then
	$show "$old_striplib $oldlib"
	$run eval "$old_striplib $oldlib" || exit $?
      fi

      # Do each command in the postinstall commands.
      eval cmds=\"$old_postinstall_cmds\"
      save_ifs="$IFS"; IFS='~'
      for cmd in $cmds; do
	IFS="$save_ifs"
	$show "$cmd"
	$run eval "$cmd" || exit $?
      done
      IFS="$save_ifs"
    done

    if test -n "$future_libdirs"; then
      $echo "$modename: warning: remember to run \`$progname --finish$future_libdirs'" 1>&2
    fi

    if test -n "$current_libdirs"; then
      # Maybe just do a dry run.
      test -n "$run" && current_libdirs=" -n$current_libdirs"
      exec_cmd='$SHELL $0 --finish$current_libdirs'
    else
      exit 0
    fi
    ;;

  # libtool finish mode
  finish)
    modename="$modename: finish"
    libdirs="$nonopt"
    admincmds=

    if test -n "$finish_cmds$finish_eval" && test -n "$libdirs"; then
      for dir
      do
	libdirs="$libdirs $dir"
      done

      for libdir in $libdirs; do
	if test -n "$finish_cmds"; then
	  # Do each command in the finish commands.
	  eval cmds=\"$finish_cmds\"
	  save_ifs="$IFS"; IFS='~'
	  for cmd in $cmds; do
	    IFS="$save_ifs"
	    $show "$cmd"
	    $run eval "$cmd" || admincmds="$admincmds
       $cmd"
	  done
	  IFS="$save_ifs"
	fi
	if test -n "$finish_eval"; then
	  # Do the single finish_eval.
	  eval cmds=\"$finish_eval\"
	  $run eval "$cmds" || admincmds="$admincmds
       $cmds"
	fi
      done
    fi

    # Exit here if they wanted silent mode.
    test "$show" = ":" && exit 0

    echo "----------------------------------------------------------------------"
    echo "Libraries have been installed in:"
    for libdir in $libdirs; do
      echo "   $libdir"
    done
    echo
    echo "If you ever happen to want to link against installed libraries"
    echo "in a given directory, LIBDIR, you must either use libtool, and"
    echo "specify the full pathname of the library, or use the \`-LLIBDIR'"
    echo "flag during linking and do at least one of the following:"
    if test -n "$shlibpath_var"; then
      echo "   - add LIBDIR to the \`$shlibpath_var' environment variable"
      echo "     during execution"
    fi
    if test -n "$runpath_var"; then
      echo "   - add LIBDIR to the \`$runpath_var' environment variable"
      echo "     during linking"
    fi
    if test -n "$hardcode_libdir_flag_spec"; then
      libdir=LIBDIR
      eval flag=\"$hardcode_libdir_flag_spec\"

      echo "   - use the \`$flag' linker flag"
    fi
    if test -n "$admincmds"; then
      echo "   - have your system administrator run these commands:$admincmds"
    fi
    if test -f /etc/ld.so.conf; then
      echo "   - have your system administrator add LIBDIR to \`/etc/ld.so.conf'"
    fi
    echo
    echo "See any operating system documentation about shared libraries for"
    echo "more information, such as the ld(1) and ld.so(8) manual pages."
    echo "----------------------------------------------------------------------"
    exit 0
    ;;

  # libtool execute mode
  execute)
    modename="$modename: execute"

    # The first argument is the command name.
    cmd="$nonopt"
    if test -z "$cmd"; then
      $echo "$modename: you must specify a COMMAND" 1>&2
      $echo "$help"
      exit 1
    fi

    # Handle -dlopen flags immediately.
    for file in $execute_dlfiles; do
      if test ! -f "$file"; then
	$echo "$modename: \`$file' is not a file" 1>&2
	$echo "$help" 1>&2
	exit 1
      fi

      dir=
      case $file in
      *.la)
	# Check to see that this really is a libtool archive.
	if (${SED} -e '2q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then :
	else
	  $echo "$modename: \`$lib' is not a valid libtool archive" 1>&2
	  $echo "$help" 1>&2
	  exit 1
	fi

	# Read the libtool library.
	dlname=
	library_names=

	# If there is no directory component, then add one.
	case $file in
	*/* | *\\*) . $file ;;
	*) . ./$file ;;
	esac

	# Skip this library if it cannot be dlopened.
	if test -z "$dlname"; then
	  # Warn if it was a shared library.
	  test -n "$library_names" && $echo "$modename: warning: \`$file' was not linked with \`-export-dynamic'"
	  continue
	fi

	dir=`$echo "X$file" | $Xsed -e 's%/[^/]*$%%'`
	test "X$dir" = "X$file" && dir=.

	if test -f "$dir/$objdir/$dlname"; then
	  dir="$dir/$objdir"
	else
	  $echo "$modename: cannot find \`$dlname' in \`$dir' or \`$dir/$objdir'" 1>&2
	  exit 1
	fi
	;;

      *.lo)
	# Just add the directory containing the .lo file.
	dir=`$echo "X$file" | $Xsed -e 's%/[^/]*$%%'`
	test "X$dir" = "X$file" && dir=.
	;;

      *)
	$echo "$modename: warning \`-dlopen' is ignored for non-libtool libraries and objects" 1>&2
	continue
	;;
      esac

      # Get the absolute pathname.
      absdir=`cd "$dir" && pwd`
      test -n "$absdir" && dir="$absdir"

      # Now add the directory to shlibpath_var.
      if eval "test -z \"\$$shlibpath_var\""; then
	eval "$shlibpath_var=\"\$dir\""
      else
	eval "$shlibpath_var=\"\$dir:\$$shlibpath_var\""
      fi
    done

    # This variable tells wrapper scripts just to set shlibpath_var
    # rather than running their programs.
    libtool_execute_magic="$magic"

    # Check if any of the arguments is a wrapper script.
    args=
    for file
    do
      case $file in
      -*) ;;
      *)
	# Do a test to see if this is really a libtool program.
	if (${SED} -e '4q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then
	  # If there is no directory component, then add one.
	  case $file in
	  */* | *\\*) . $file ;;
	  *) . ./$file ;;
	  esac

	  # Transform arg to wrapped name.
	  file="$progdir/$program"
	fi
	;;
      esac
      # Quote arguments (to preserve shell metacharacters).
      file=`$echo "X$file" | $Xsed -e "$sed_quote_subst"`
      args="$args \"$file\""
    done

    if test -z "$run"; then
      if test -n "$shlibpath_var"; then
	# Export the shlibpath_var.
	eval "export $shlibpath_var"
      fi

      # Restore saved enviroment variables
      if test "${save_LC_ALL+set}" = set; then
	LC_ALL="$save_LC_ALL"; export LC_ALL
      fi
      if test "${save_LANG+set}" = set; then
	LANG="$save_LANG"; export LANG
      fi

      # Now prepare to actually exec the command.
      exec_cmd="\$cmd$args"
    else
      # Display what would be done.
      if test -n "$shlibpath_var"; then
	eval "\$echo \"\$shlibpath_var=\$$shlibpath_var\""
	$echo "export $shlibpath_var"
      fi
      $echo "$cmd$args"
      exit 0
    fi
    ;;

  # libtool clean and uninstall mode
  clean | uninstall)
    modename="$modename: $mode"
    rm="$nonopt"
    files=
    rmforce=
    exit_status=0

    # This variable tells wrapper scripts just to set variables rather
    # than running their programs.
    libtool_install_magic="$magic"

    for arg
    do
      case $arg in
      -f) rm="$rm $arg"; rmforce=yes ;;
      -*) rm="$rm $arg" ;;
      *) files="$files $arg" ;;
      esac
    done

    if test -z "$rm"; then
      $echo "$modename: you must specify an RM program" 1>&2
      $echo "$help" 1>&2
      exit 1
    fi

    rmdirs=

    for file in $files; do
      dir=`$echo "X$file" | $Xsed -e 's%/[^/]*$%%'`
      if test "X$dir" = "X$file"; then
	dir=.
	objdir="$objdir"
      else
	objdir="$dir/$objdir"
      fi
      name=`$echo "X$file" | $Xsed -e 's%^.*/%%'`
      test $mode = uninstall && objdir="$dir"

      # Remember objdir for removal later, being careful to avoid duplicates
      if test $mode = clean; then
	case " $rmdirs " in
	  *" $objdir "*) ;;
	  *) rmdirs="$rmdirs $objdir" ;;
	esac
      fi

      # Don't error if the file doesn't exist and rm -f was used.
      if (test -L "$file") >/dev/null 2>&1 \
	|| (test -h "$file") >/dev/null 2>&1 \
	|| test -f "$file"; then
	:
      elif test -d "$file"; then
	exit_status=1
	continue
      elif test "$rmforce" = yes; then
	continue
      fi

      rmfiles="$file"

      case $name in
      *.la)
	# Possibly a libtool archive, so verify it.
	if (${SED} -e '2q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then
	  . $dir/$name

	  # Delete the libtool libraries and symlinks.
	  for n in $library_names; do
	    rmfiles="$rmfiles $objdir/$n"
	  done
	  test -n "$old_library" && rmfiles="$rmfiles $objdir/$old_library"
	  test $mode = clean && rmfiles="$rmfiles $objdir/$name $objdir/${name}i"

	  if test $mode = uninstall; then
	    if test -n "$library_names"; then
	      # Do each command in the postuninstall commands.
	      eval cmds=\"$postuninstall_cmds\"
	      save_ifs="$IFS"; IFS='~'
	      for cmd in $cmds; do
		IFS="$save_ifs"
		$show "$cmd"
		$run eval "$cmd"
		if test $? != 0 && test "$rmforce" != yes; then
		  exit_status=1
		fi
	      done
	      IFS="$save_ifs"
	    fi

	    if test -n "$old_library"; then
	      # Do each command in the old_postuninstall commands.
	      eval cmds=\"$old_postuninstall_cmds\"
	      save_ifs="$IFS"; IFS='~'
	      for cmd in $cmds; do
		IFS="$save_ifs"
		$show "$cmd"
		$run eval "$cmd"
		if test $? != 0 && test "$rmforce" != yes; then
		  exit_status=1
		fi
	      done
	      IFS="$save_ifs"
	    fi
	    # FIXME: should reinstall the best remaining shared library.
	  fi
	fi
	;;

      *.lo)
	# Possibly a libtool object, so verify it.
	if (sed -e '2q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then

	  # Read the .lo file
	  . $dir/$name

	  # Add PIC object to the list of files to remove.
	  if test -n "$pic_object" \
	     && test "$pic_object" != none; then
	    rmfiles="$rmfiles $dir/$pic_object"
	  fi

	  # Add non-PIC object to the list of files to remove.
	  if test -n "$non_pic_object" \
	     && test "$non_pic_object" != none; then
	    rmfiles="$rmfiles $dir/$non_pic_object"
	  fi
	fi
	;;

      *)
	# Do a test to see if this is a libtool program.
	if test $mode = clean &&
	   (${SED} -e '4q' $file | egrep "^# Generated by .*$PACKAGE") >/dev/null 2>&1; then
	  relink_command=
	  . $dir/$file

	  rmfiles="$rmfiles $objdir/$name $objdir/${name}S.${objext}"
	  if test "$fast_install" = yes && test -n "$relink_command"; then
	    rmfiles="$rmfiles $objdir/lt-$name"
	  fi
	fi
	;;
      esac
      $show "$rm $rmfiles"
      $run $rm $rmfiles || exit_status=1
    done

    # Try to remove the ${objdir}s in the directories where we deleted files
    for dir in $rmdirs; do
      if test -d "$dir"; then
	$show "rmdir $dir"
	$run rmdir $dir >/dev/null 2>&1
      fi
    done

    exit $exit_status
    ;;

  "")
    $echo "$modename: you must specify a MODE" 1>&2
    $echo "$generic_help" 1>&2
    exit 1
    ;;
  esac

  if test -z "$exec_cmd"; then
    $echo "$modename: invalid operation mode \`$mode'" 1>&2
    $echo "$generic_help" 1>&2
    exit 1
  fi
fi # test -z "$show_help"

if test -n "$exec_cmd"; then
  eval exec $exec_cmd
  exit 1
fi

# We need to display help for each of the modes.
case $mode in
"") $echo \
"Usage: $modename [OPTION]... [MODE-ARG]...

Provide generalized library-building support services.

    --config          show all configuration variables
    --debug           enable verbose shell tracing
-n, --dry-run         display commands without modifying any files
    --features        display basic configuration information and exit
    --finish          same as \`--mode=finish'
    --help            display this help message and exit
    --mode=MODE       use operation mode MODE [default=inferred from MODE-ARGS]
    --quiet           same as \`--silent'
    --silent          don't print informational messages
    --version         print version information

MODE must be one of the following:

      clean           remove files from the build directory
      compile         compile a source file into a libtool object
      execute         automatically set library path, then run a program
      finish          complete the installation of libtool libraries
      install         install libraries or executables
      link            create a library or an executable
      uninstall       remove libraries from an installed directory

MODE-ARGS vary depending on the MODE.  Try \`$modename --help --mode=MODE' for
a more detailed description of MODE."
  exit 0
  ;;

clean)
  $echo \
"Usage: $modename [OPTION]... --mode=clean RM [RM-OPTION]... FILE...

Remove files from the build directory.

RM is the name of the program to use to delete files associated with each FILE
(typically \`/bin/rm').  RM-OPTIONS are options (such as \`-f') to be passed
to RM.

If FILE is a libtool library, object or program, all the files associated
with it are deleted. Otherwise, only FILE itself is deleted using RM."
  ;;

compile)
  $echo \
"Usage: $modename [OPTION]... --mode=compile COMPILE-COMMAND... SOURCEFILE

Compile a source file into a libtool library object.

This mode accepts the following additional options:

  -o OUTPUT-FILE    set the output file name to OUTPUT-FILE
  -prefer-pic       try to building PIC objects only
  -prefer-non-pic   try to building non-PIC objects only
  -static           always build a \`.o' file suitable for static linking

COMPILE-COMMAND is a command to be used in creating a \`standard' object file
from the given SOURCEFILE.

The output file name is determined by removing the directory component from
SOURCEFILE, then substituting the C source code suffix \`.c' with the
library object suffix, \`.lo'."
  ;;

execute)
  $echo \
"Usage: $modename [OPTION]... --mode=execute COMMAND [ARGS]...

Automatically set library path, then run a program.

This mode accepts the following additional options:

  -dlopen FILE      add the directory containing FILE to the library path

This mode sets the library path environment variable according to \`-dlopen'
flags.

If any of the ARGS are libtool executable wrappers, then they are translated
into their corresponding uninstalled binary, and any of their required library
directories are added to the library path.

Then, COMMAND is executed, with ARGS as arguments."
  ;;

finish)
  $echo \
"Usage: $modename [OPTION]... --mode=finish [LIBDIR]...

Complete the installation of libtool libraries.

Each LIBDIR is a directory that contains libtool libraries.

The commands that this mode executes may require superuser privileges.  Use
the \`--dry-run' option if you just want to see what would be executed."
  ;;

install)
  $echo \
"Usage: $modename [OPTION]... --mode=install INSTALL-COMMAND...

Install executables or libraries.

INSTALL-COMMAND is the installation command.  The first component should be
either the \`install' or \`cp' program.

The rest of the components are interpreted as arguments to that command (only
BSD-compatible install options are recognized)."
  ;;

link)
  $echo \
"Usage: $modename [OPTION]... --mode=link LINK-COMMAND...

Link object files or libraries together to form another library, or to
create an executable program.

LINK-COMMAND is a command using the C compiler that you would use to create
a program from several object files.

The following components of LINK-COMMAND are treated specially:

  -all-static       do not do any dynamic linking at all
  -avoid-version    do not add a version suffix if possible
  -dlopen FILE      \`-dlpreopen' FILE if it cannot be dlopened at runtime
  -dlpreopen FILE   link in FILE and add its symbols to lt_preloaded_symbols
  -export-dynamic   allow symbols from OUTPUT-FILE to be resolved with dlsym(3)
  -export-symbols SYMFILE
		    try to export only the symbols listed in SYMFILE
  -export-symbols-regex REGEX
		    try to export only the symbols matching REGEX
  -LLIBDIR          search LIBDIR for required installed libraries
  -lNAME            OUTPUT-FILE requires the installed library libNAME
  -module           build a library that can dlopened
  -no-fast-install  disable the fast-install mode
  -no-install       link a not-installable executable
  -no-undefined     declare that a library does not refer to external symbols
  -o OUTPUT-FILE    create OUTPUT-FILE from the specified objects
  -release RELEASE  specify package release information
  -rpath LIBDIR     the created library will eventually be installed in LIBDIR
  -R[ ]LIBDIR       add LIBDIR to the runtime path of programs and libraries
  -static           do not do any dynamic linking of libtool libraries
  -version-info CURRENT[:REVISION[:AGE]]
		    specify library version info [each variable defaults to 0]

All other options (arguments beginning with \`-') are ignored.

Every other argument is treated as a filename.  Files ending in \`.la' are
treated as uninstalled libtool libraries, other files are standard or library
object files.

If the OUTPUT-FILE ends in \`.la', then a libtool library is created,
only library objects (\`.lo' files) may be specified, and \`-rpath' is
required, except when creating a convenience library.

If OUTPUT-FILE ends in \`.a' or \`.lib', then a standard library is created
using \`ar' and \`ranlib', or on Windows using \`lib'.

If OUTPUT-FILE ends in \`.lo' or \`.${objext}', then a reloadable object file
is created, otherwise an executable program is created."
  ;;

uninstall)
  $echo \
"Usage: $modename [OPTION]... --mode=uninstall RM [RM-OPTION]... FILE...

Remove libraries from an installation directory.

RM is the name of the program to use to delete files associated with each FILE
(typically \`/bin/rm').  RM-OPTIONS are options (such as \`-f') to be passed
to RM.

If FILE is a libtool library, all the files associated with it are deleted.
Otherwise, only FILE itself is deleted using RM."
  ;;

*)
  $echo "$modename: invalid operation mode \`$mode'" 1>&2
  $echo "$help" 1>&2
  exit 1
  ;;
esac

echo
$echo "Try \`$modename --help' for more information about other modes."

exit 0

# Local Variables:
# mode:shell-script
# sh-indentation:2
# End:
