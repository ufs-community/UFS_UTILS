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

if [[ -d /lfs5 ]] ; then
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
elif [[ -d /opt/spack-stack ]] ; then
    # We are using a container 
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        echo load the module command 1>&2
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    target=container
    module purge
elif [[ -d /scratch3 ]] ; then
    # We are on NOAA Ursa
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        echo load the module command 1>&2
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    target=ursa
    module purge
elif [[ "$(hostname)" == "gaea6"* || "$(hostname)" =~ c6n[0-9]+ ]] && [[ -d /gpfs/f6 ]] ; then
    target=gaeac6
    source /opt/cray/pe/lmod/8.7.31/init/$__ms_shell
elif [[ "$(hostname)" =~ "Orion" || "$(hostname)" =~ "orion" ]]; then
    target="orion"
    module purge
elif [[ "$(hostname)" =~ "hercules" || "$(hostname)" =~ "Hercules" ]]; then
    target="hercules"
    module purge
elif [[ -d /work/00315 && -d /scratch/00315 ]] ; then
    target=stampede
    module purge
else
    if [[ ! -v PW_CSP ]]; then
        set +x
        echo FATAL ERROR: UNKNOWN PLATFORM 1>&2; exit 99
    elif [[ -z "${PW_CSP}" ]]; then
        set +x
        echo FATAL ERROR: UNKNOWN PLATFORM 1>&2; exit 99
    else
        if [[ "${PW_CSP}" == "aws" || "${PW_CSP}" == "azure" || "${PW_CSP}" == "google" ]]; then
            target=noaacloud
            module purge
        else
            set +x
            echo FATAL ERROR: UNKNOWN PLATFORM 1>&2; exit 99
        fi
    fi
fi

unset __ms_shell
unset __ms_ksh_test
unset __ms_bash_test
unset $__ms_function_name
unset __ms_function_name
