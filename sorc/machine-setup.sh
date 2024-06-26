# Create a test function for sh vs. bash detection.  The name is
# randomly generated to reduce the chances of name collision.
__ms_function_name="setup__test_function__$$"
eval "$__ms_function_name() { /bin/true ; }"

# Determine which shell we are using
__ms_ksh_test=$( eval '__text="text" ; if [[ $__text =~ ^(t).* ]] ; then printf "%s" ${.sh.match[1]} ; fi' 2> /dev/null | cat )
__ms_bash_test=$( eval 'if ( set | grep '$__ms_function_name' | grep -v name > /dev/null 2>&1 ) ; then echo t ; fi ' 2> /dev/null | cat )

if [[ ! -z "$__ms_ksh_test" ]] ; then
    __ms_shell=ksh
elif [[ ! -z "$__ms_bash_test" ]] ; then
    __ms_shell=bash
else
    # Not bash or ksh, so assume sh.
    __ms_shell=sh
fi

target=""
USERNAME=`echo $LOGNAME | awk '{ print tolower($0)'}`

if [[ -d /lfs4 ]] ; then
    # We are on NOAA Jet
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        echo load the module command 1>&2
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    target=jet
    module purge
elif [[ -d /lfs/h1 ]] ; then
    target=wcoss2
    module reset
elif [[ -d /scratch1 ]] ; then
    # We are on NOAA Hera
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        echo load the module command 1>&2
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    target=hera
    module purge
elif [[ -d /gpfs && -d /ncrc ]] ; then
    # We are on GAEA.
    if ( ! eval module help > /dev/null 2>&1 ) ; then
      # We cannot simply load the module command.  The GAEA
      # /etc/profile modifies a number of module-related variables
      # before loading the module command.  Without those variables,
      # the module command fails.  Hence we actually have to source
      # /etc/profile here.
      source /etc/profile
    fi
    module reset
    target=gaea
elif [[ "$(hostname)" =~ "Orion" || "$(hostname)" =~ "orion" ]]; then
    target="orion"
    module purge
elif [[ "$(hostname)" =~ "hercules" || "$(hostname)" =~ "Hercules" ]]; then
    target="hercules"
    module purge
elif [[ -d /work/00315 && -d /scratch/00315 ]] ; then
    target=stampede
    module purge
elif [[ -d /data/prod ]] ; then
    # We are on SSEC S4
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        echo load the module command 1>&2
        source /usr/share/lmod/lmod/init/$__ms_shell
    fi
    target=s4
    module purge
else
    if [[ ! -v PW_CSP ]]; then
        echo WARNING: UNKNOWN PLATFORM 1>&2; exit 99
    elif [[ -z "${PW_CSP}" ]]; then
        echo WARNING: UNKNOWN PLATFORM 1>&2; exit 99
    else
        if [[ "${PW_CSP}" == "aws" || "${PW_CSP}" == "azure" || "${PW_CSP}" == "google" ]]; then
            target=noaacloud
            module purge
        else
            echo WARNING: UNKNOWN PLATFORM 1>&2; exit 99
        fi
    fi
fi

unset __ms_shell
unset __ms_ksh_test
unset __ms_bash_test
unset $__ms_function_name
unset __ms_function_name
