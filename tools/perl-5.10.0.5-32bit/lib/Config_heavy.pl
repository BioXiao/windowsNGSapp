# This file was created by configpm when Perl was built. Any changes
# made to this file will be lost the next time perl is built.

package Config;
use strict;
# use warnings; Pulls in Carp
# use vars pulls in Carp
### Configured by: Adam.Kennedy@SYD1D28264.ce.corp
### Target system: WIN32 

our $summary = <<'!END!';
Summary of my $package (revision $revision $version_patchlevel_string) configuration:
  Platform:
    osname=$osname, osvers=$osvers, archname=$archname
    uname='$myuname'
    config_args='$config_args'
    hint=$hint, useposix=$useposix, d_sigaction=$d_sigaction
    useithreads=$useithreads, usemultiplicity=$usemultiplicity
    useperlio=$useperlio, d_sfio=$d_sfio, uselargefiles=$uselargefiles, usesocks=$usesocks
    use64bitint=$use64bitint, use64bitall=$use64bitall, uselongdouble=$uselongdouble
    usemymalloc=$usemymalloc, bincompat5005=undef
  Compiler:
    cc='$cc', ccflags ='$ccflags',
    optimize='$optimize',
    cppflags='$cppflags'
    ccversion='$ccversion', gccversion='$gccversion', gccosandvers='$gccosandvers'
    intsize=$intsize, longsize=$longsize, ptrsize=$ptrsize, doublesize=$doublesize, byteorder=$byteorder
    d_longlong=$d_longlong, longlongsize=$longlongsize, d_longdbl=$d_longdbl, longdblsize=$longdblsize
    ivtype='$ivtype', ivsize=$ivsize, nvtype='$nvtype', nvsize=$nvsize, Off_t='$lseektype', lseeksize=$lseeksize
    alignbytes=$alignbytes, prototype=$prototype
  Linker and Libraries:
    ld='$ld', ldflags ='$ldflags'
    libpth=$libpth
    libs=$libs
    perllibs=$perllibs
    libc=$libc, so=$so, useshrplib=$useshrplib, libperl=$libperl
    gnulibc_version='$gnulibc_version'
  Dynamic Linking:
    dlsrc=$dlsrc, dlext=$dlext, d_dlsymun=$d_dlsymun, ccdlflags='$ccdlflags'
    cccdlflags='$cccdlflags', lddlflags='$lddlflags'

!END!
my $summary_expanded;

sub myconfig {
    return $summary_expanded if $summary_expanded;
    ($summary_expanded = $summary) =~ s{\$(\w+)}
		 { my $c = $Config::Config{$1}; defined($c) ? $c : 'undef' }ge;
    $summary_expanded;
}

local *_ = \my $a;
$_ = <<'!END!';
Author=''
CONFIG='true'
Date='$Date'
Header=''
Id='$Id'
Locker=''
Log='$Log'
Mcc='Mcc'
PATCHLEVEL='10'
PERL_API_REVISION='5'
PERL_API_SUBVERSION='0'
PERL_API_VERSION='10'
PERL_CONFIG_SH='true'
PERL_PATCHLEVEL=''
PERL_REVISION='5'
PERL_SUBVERSION='0'
PERL_VERSION='10'
RCSfile='$RCSfile'
Revision='$Revision'
SUBVERSION='0'
Source=''
State=''
_a='.a'
_exe='.exe'
_o='.o'
afs='false'
afsroot='/afs'
alignbytes='8'
ansi2knr=''
aphostname=''
api_revision='5'
api_subversion='0'
api_version='10'
api_versionstring='5.10.0'
ar='ar'
archlib='C:\strawberry\perl\lib'
archlibexp='C:\strawberry\perl\lib'
archname64=''
archname='MSWin32-x86-multi-thread'
archobjs=''
asctime_r_proto='0'
awk='awk'
baserev='5'
bash=''
bin='C:\strawberry\perl\bin'
binexp='C:\strawberry\perl\bin'
bison=''
byacc='byacc'
byteorder='1234'
c=''
castflags='0'
cat='type'
cc='gcc'
cccdlflags=' '
ccdlflags=' '
ccflags=' -s -O2 -DWIN32 -DHAVE_DES_FCRYPT  -DPERL_IMPLICIT_CONTEXT -DPERL_IMPLICIT_SYS -fno-strict-aliasing -DPERL_MSVCRT_READFIX'
ccflags_uselargefiles=''
ccname='gcc'
ccsymbols=''
ccversion=''
cf_by='Adam.Kennedy'
cf_email='Adam.Kennedy@SYD1D28264.ce.corp'
cf_time='Thu Apr 30 15:18:53 2009'
chgrp=''
chmod=''
chown=''
clocktype='clock_t'
comm=''
compress=''
contains='grep'
cp='copy'
cpio=''
cpp='gcc -E'
cpp_stuff='42'
cppccsymbols=''
cppflags='-DWIN32'
cpplast=''
cppminus='-'
cpprun='gcc -E'
cppstdin='gcc -E'
cppsymbols=''
crypt_r_proto='0'
cryptlib=''
csh='undef'
ctermid_r_proto='0'
ctime_r_proto='0'
d_Gconvert='sprintf((b),"%.*g",(n),(x))'
d_PRIEUldbl='undef'
d_PRIFUldbl='undef'
d_PRIGUldbl='undef'
d_PRIXU64='undef'
d_PRId64='undef'
d_PRIeldbl='undef'
d_PRIfldbl='undef'
d_PRIgldbl='undef'
d_PRIi64='undef'
d_PRIo64='undef'
d_PRIu64='undef'
d_PRIx64='undef'
d_SCNfldbl='undef'
d__fwalk='undef'
d_access='define'
d_accessx='undef'
d_aintl='undef'
d_alarm='define'
d_archlib='define'
d_asctime_r='undef'
d_atolf='undef'
d_atoll='undef'
d_attribute_format='undef'
d_attribute_malloc='undef'
d_attribute_nonnull='undef'
d_attribute_noreturn='undef'
d_attribute_pure='undef'
d_attribute_unused='undef'
d_attribute_warn_unused_result='undef'
d_bcmp='undef'
d_bcopy='undef'
d_bsd='define'
d_bsdgetpgrp='undef'
d_bsdsetpgrp='undef'
d_builtin_choose_expr='undef'
d_builtin_expect='undef'
d_bzero='undef'
d_c99_variadic_macros='undef'
d_casti32='define'
d_castneg='define'
d_charvspr='undef'
d_chown='undef'
d_chroot='undef'
d_chsize='define'
d_class='undef'
d_clearenv='undef'
d_closedir='define'
d_cmsghdr_s='undef'
d_const='define'
d_copysignl='undef'
d_cplusplus='undef'
d_crypt='define'
d_crypt_r='undef'
d_csh='undef'
d_ctermid='undef'
d_ctermid_r='undef'
d_ctime_r='undef'
d_cuserid='undef'
d_dbl_dig='define'
d_dbminitproto='undef'
d_difftime='define'
d_dir_dd_fd='undef'
d_dirfd='undef'
d_dirnamlen='define'
d_dlerror='define'
d_dlopen='define'
d_dlsymun='undef'
d_dosuid='undef'
d_drand48_r='undef'
d_drand48proto='undef'
d_dup2='define'
d_eaccess='undef'
d_endgrent='undef'
d_endgrent_r='undef'
d_endhent='undef'
d_endhostent_r='undef'
d_endnent='undef'
d_endnetent_r='undef'
d_endpent='undef'
d_endprotoent_r='undef'
d_endpwent='undef'
d_endpwent_r='undef'
d_endsent='undef'
d_endservent_r='undef'
d_eofnblk='define'
d_eunice='undef'
d_faststdio='define'
d_fchdir='undef'
d_fchmod='undef'
d_fchown='undef'
d_fcntl='undef'
d_fcntl_can_lock='undef'
d_fd_macros='define'
d_fd_set='define'
d_fds_bits='define'
d_fgetpos='define'
d_finite='undef'
d_finitel='undef'
d_flexfnam='define'
d_flock='define'
d_flockproto='define'
d_fork='undef'
d_fp_class='undef'
d_fpathconf='undef'
d_fpclass='undef'
d_fpclassify='undef'
d_fpclassl='undef'
d_fpos64_t='undef'
d_frexpl='undef'
d_fs_data_s='undef'
d_fseeko='undef'
d_fsetpos='define'
d_fstatfs='undef'
d_fstatvfs='undef'
d_fsync='undef'
d_ftello='undef'
d_ftime='define'
d_futimes='undef'
d_getcwd='define'
d_getespwnam='undef'
d_getfsstat='undef'
d_getgrent='undef'
d_getgrent_r='undef'
d_getgrgid_r='undef'
d_getgrnam_r='undef'
d_getgrps='undef'
d_gethbyaddr='define'
d_gethbyname='define'
d_gethent='undef'
d_gethname='define'
d_gethostbyaddr_r='undef'
d_gethostbyname_r='undef'
d_gethostent_r='undef'
d_gethostprotos='define'
d_getitimer='undef'
d_getlogin='define'
d_getlogin_r='undef'
d_getmnt='undef'
d_getmntent='undef'
d_getnbyaddr='undef'
d_getnbyname='undef'
d_getnent='undef'
d_getnetbyaddr_r='undef'
d_getnetbyname_r='undef'
d_getnetent_r='undef'
d_getnetprotos='undef'
d_getpagsz='undef'
d_getpbyname='define'
d_getpbynumber='define'
d_getpent='undef'
d_getpgid='undef'
d_getpgrp2='undef'
d_getpgrp='undef'
d_getppid='undef'
d_getprior='undef'
d_getprotobyname_r='undef'
d_getprotobynumber_r='undef'
d_getprotoent_r='undef'
d_getprotoprotos='define'
d_getprpwnam='undef'
d_getpwent='undef'
d_getpwent_r='undef'
d_getpwnam_r='undef'
d_getpwuid_r='undef'
d_getsbyname='define'
d_getsbyport='define'
d_getsent='undef'
d_getservbyname_r='undef'
d_getservbyport_r='undef'
d_getservent_r='undef'
d_getservprotos='define'
d_getspnam='undef'
d_getspnam_r='undef'
d_gettimeod='define'
d_gmtime_r='undef'
d_gnulibc='undef'
d_grpasswd='undef'
d_hasmntopt='undef'
d_htonl='define'
d_ilogbl='undef'
d_inc_version_list='undef'
d_index='undef'
d_inetaton='undef'
d_int64_t='undef'
d_isascii='define'
d_isfinite='undef'
d_isinf='undef'
d_isnan='define'
d_isnanl='undef'
d_killpg='define'
d_lchown='undef'
d_ldbl_dig='define'
d_libm_lib_version='undef'
d_link='define'
d_localtime_r='undef'
d_localtime_r_needs_tzset='undef'
d_locconv='define'
d_lockf='undef'
d_longdbl='define'
d_longlong='undef'
d_lseekproto='define'
d_lstat='undef'
d_madvise='undef'
d_malloc_good_size='undef'
d_malloc_size='undef'
d_mblen='define'
d_mbstowcs='define'
d_mbtowc='define'
d_memchr='define'
d_memcmp='define'
d_memcpy='define'
d_memmove='define'
d_memset='define'
d_mkdir='define'
d_mkdtemp='undef'
d_mkfifo='undef'
d_mkstemp='undef'
d_mkstemps='undef'
d_mktime='define'
d_mmap='undef'
d_modfl='undef'
d_modfl_pow32_bug='undef'
d_modflproto='undef'
d_mprotect='undef'
d_msg='undef'
d_msg_ctrunc='undef'
d_msg_dontroute='undef'
d_msg_oob='undef'
d_msg_peek='undef'
d_msg_proxy='undef'
d_msgctl='undef'
d_msgget='undef'
d_msghdr_s='undef'
d_msgrcv='undef'
d_msgsnd='undef'
d_msync='undef'
d_munmap='undef'
d_mymalloc='undef'
d_nice='undef'
d_nl_langinfo='undef'
d_nv_preserves_uv='define'
d_nv_zero_is_allbits_zero='define'
d_off64_t='undef'
d_old_pthread_create_joinable='undef'
d_oldpthreads='undef'
d_oldsock='undef'
d_open3='undef'
d_pathconf='undef'
d_pause='define'
d_perl_otherlibdirs='undef'
d_phostname='undef'
d_pipe='define'
d_poll='undef'
d_portable='define'
d_printf_format_null='undef'
d_procselfexe='undef'
d_pseudofork='define'
d_pthread_atfork='undef'
d_pthread_attr_setscope='undef'
d_pthread_yield='undef'
d_pwage='undef'
d_pwchange='undef'
d_pwclass='undef'
d_pwcomment='undef'
d_pwexpire='undef'
d_pwgecos='undef'
d_pwpasswd='undef'
d_pwquota='undef'
d_qgcvt='undef'
d_quad='undef'
d_random_r='undef'
d_readdir64_r='undef'
d_readdir='define'
d_readdir_r='undef'
d_readlink='undef'
d_readv='undef'
d_recvmsg='undef'
d_rename='define'
d_rewinddir='define'
d_rmdir='define'
d_safebcpy='undef'
d_safemcpy='undef'
d_sanemcmp='define'
d_sbrkproto='undef'
d_scalbnl='undef'
d_sched_yield='undef'
d_scm_rights='undef'
d_seekdir='define'
d_select='define'
d_sem='undef'
d_semctl='undef'
d_semctl_semid_ds='undef'
d_semctl_semun='undef'
d_semget='undef'
d_semop='undef'
d_sendmsg='undef'
d_setegid='undef'
d_seteuid='undef'
d_setgrent='undef'
d_setgrent_r='undef'
d_setgrps='undef'
d_sethent='undef'
d_sethostent_r='undef'
d_setitimer='undef'
d_setlinebuf='undef'
d_setlocale='define'
d_setlocale_r='undef'
d_setnent='undef'
d_setnetent_r='undef'
d_setpent='undef'
d_setpgid='undef'
d_setpgrp2='undef'
d_setpgrp='undef'
d_setprior='undef'
d_setproctitle='undef'
d_setprotoent_r='undef'
d_setpwent='undef'
d_setpwent_r='undef'
d_setregid='undef'
d_setresgid='undef'
d_setresuid='undef'
d_setreuid='undef'
d_setrgid='undef'
d_setruid='undef'
d_setsent='undef'
d_setservent_r='undef'
d_setsid='undef'
d_setvbuf='define'
d_sfio='undef'
d_shm='undef'
d_shmat='undef'
d_shmatprototype='undef'
d_shmctl='undef'
d_shmdt='undef'
d_shmget='undef'
d_sigaction='undef'
d_signbit='undef'
d_sigprocmask='undef'
d_sigsetjmp='undef'
d_sitearch='define'
d_snprintf='define'
d_sockatmark='undef'
d_sockatmarkproto='undef'
d_socket='define'
d_socklen_t='undef'
d_sockpair='undef'
d_socks5_init='undef'
d_sprintf_returns_strlen='define'
d_sqrtl='undef'
d_srand48_r='undef'
d_srandom_r='undef'
d_sresgproto='undef'
d_sresuproto='undef'
d_statblks='undef'
d_statfs_f_flags='undef'
d_statfs_s='undef'
d_statvfs='undef'
d_stdio_cnt_lval='define'
d_stdio_ptr_lval='define'
d_stdio_ptr_lval_nochange_cnt='define'
d_stdio_ptr_lval_sets_cnt='undef'
d_stdio_stream_array='undef'
d_stdiobase='define'
d_stdstdio='define'
d_strchr='define'
d_strcoll='define'
d_strctcpy='define'
d_strerrm='strerror(e)'
d_strerror='define'
d_strerror_r='undef'
d_strftime='define'
d_strlcat='undef'
d_strlcpy='undef'
d_strtod='define'
d_strtol='define'
d_strtold='undef'
d_strtoll='undef'
d_strtoq='undef'
d_strtoul='define'
d_strtoull='undef'
d_strtouq='undef'
d_strxfrm='define'
d_suidsafe='undef'
d_symlink='undef'
d_syscall='undef'
d_syscallproto='undef'
d_sysconf='undef'
d_sysernlst=''
d_syserrlst='define'
d_system='define'
d_tcgetpgrp='undef'
d_tcsetpgrp='undef'
d_telldir='define'
d_telldirproto='define'
d_time='define'
d_times='define'
d_tm_tm_gmtoff='undef'
d_tm_tm_zone='undef'
d_tmpnam_r='undef'
d_truncate='undef'
d_ttyname_r='undef'
d_tzname='define'
d_u32align='define'
d_ualarm='undef'
d_umask='define'
d_uname='define'
d_union_semun='define'
d_unordered='undef'
d_unsetenv='undef'
d_usleep='undef'
d_usleepproto='undef'
d_ustat='undef'
d_vendorarch='undef'
d_vendorbin='undef'
d_vendorlib='undef'
d_vendorscript='undef'
d_vfork='undef'
d_void_closedir='undef'
d_voidsig='define'
d_voidtty=''
d_volatile='define'
d_vprintf='define'
d_vsnprintf='define'
d_wait4='undef'
d_waitpid='define'
d_wcstombs='define'
d_wctomb='define'
d_writev='undef'
d_xenix='undef'
date='date'
db_hashtype='int'
db_prefixtype='int'
db_version_major='0'
db_version_minor='0'
db_version_patch='0'
defvoidused='15'
direntrytype='struct direct'
dlext='dll'
dlsrc='dl_win32.xs'
doublesize='8'
drand01='(rand()/(double)((unsigned)1<<RANDBITS))'
drand48_r_proto='0'
dynamic_ext='B Compress/Raw/Zlib Cwd Data/Dumper Devel/DProf Devel/PPPort Devel/Peek Digest/MD5 Digest/SHA Encode Fcntl File/Glob Filter/Util/Call Hash/Util Hash/Util/FieldHash IO List/Util MIME/Base64 Math/BigInt/FastCalc Opcode POSIX PerlIO/encoding PerlIO/scalar PerlIO/via SDBM_File Socket Storable Sys/Hostname Text/Soundex Time/HiRes Time/Piece Unicode/Normalize Win32 Win32API/File XS/APItest XS/Typemap attrs re threads threads/shared'
eagain='EAGAIN'
ebcdic='undef'
echo='echo'
egrep='egrep'
emacs=''
endgrent_r_proto='0'
endhostent_r_proto='0'
endnetent_r_proto='0'
endprotoent_r_proto='0'
endpwent_r_proto='0'
endservent_r_proto='0'
eunicefix=':'
exe_ext='.exe'
expr='expr'
extensions='B Compress/Raw/Zlib Compress/Zlib Cwd Data/Dumper Devel/DProf Devel/PPPort Devel/Peek Digest/MD5 Digest/SHA Encode Errno Fcntl File/Glob Filter/Util/Call Hash/Util Hash/Util/FieldHash IO IO_Compress_Base IO_Compress_Zlib List/Util MIME/Base64 Math/BigInt/FastCalc Opcode POSIX PerlIO/encoding PerlIO/scalar PerlIO/via SDBM_File Socket Storable Sys/Hostname Text/Soundex Time/HiRes Time/Piece Unicode/Normalize Win32 Win32API/File Win32CORE XS/APItest XS/Typemap attrs re threads threads/shared'
extras=''
fflushNULL='define'
fflushall='undef'
find='find'
firstmakefile='makefile'
flex=''
fpossize='8'
fpostype='fpos_t'
freetype='void'
from=':'
full_ar=''
full_csh=''
full_sed=''
gccansipedantic=''
gccosandvers=''
gccversion='3.4.5'
getgrent_r_proto='0'
getgrgid_r_proto='0'
getgrnam_r_proto='0'
gethostbyaddr_r_proto='0'
gethostbyname_r_proto='0'
gethostent_r_proto='0'
getlogin_r_proto='0'
getnetbyaddr_r_proto='0'
getnetbyname_r_proto='0'
getnetent_r_proto='0'
getprotobyname_r_proto='0'
getprotobynumber_r_proto='0'
getprotoent_r_proto='0'
getpwent_r_proto='0'
getpwnam_r_proto='0'
getpwuid_r_proto='0'
getservbyname_r_proto='0'
getservbyport_r_proto='0'
getservent_r_proto='0'
getspnam_r_proto='0'
gidformat='"ld"'
gidsign='-1'
gidsize='4'
gidtype='gid_t'
glibpth='/usr/shlib  /lib/pa1.1 /usr/lib/large /lib /usr/lib /usr/lib/386 /lib/386 /lib/large /usr/lib/small /lib/small /usr/ccs/lib /usr/ucblib /usr/shlib '
gmake='gmake'
gmtime_r_proto='0'
gnulibc_version=''
grep='grep'
groupcat=''
groupstype='gid_t'
gzip='gzip'
h_fcntl='false'
h_sysfile='true'
hint='recommended'
hostcat='ypcat hosts'
html1dir=''
html1direxp=''
html3dir=''
html3direxp=''
i16size='2'
i16type='short'
i32size='4'
i32type='long'
i64size='8'
i64type='long long'
i8size='1'
i8type='char'
i_arpainet='define'
i_bsdioctl=''
i_crypt='undef'
i_db='undef'
i_dbm='undef'
i_dirent='define'
i_dld='undef'
i_dlfcn='define'
i_fcntl='define'
i_float='define'
i_fp='undef'
i_fp_class='undef'
i_gdbm='undef'
i_grp='undef'
i_ieeefp='undef'
i_inttypes='undef'
i_langinfo='undef'
i_libutil='undef'
i_limits='define'
i_locale='define'
i_machcthr='undef'
i_malloc='define'
i_math='define'
i_memory='undef'
i_mntent='undef'
i_ndbm='undef'
i_netdb='undef'
i_neterrno='undef'
i_netinettcp='undef'
i_niin='undef'
i_poll='undef'
i_prot='undef'
i_pthread='undef'
i_pwd='undef'
i_rpcsvcdbm='define'
i_sfio='undef'
i_sgtty='undef'
i_shadow='undef'
i_socks='undef'
i_stdarg='define'
i_stddef='define'
i_stdlib='define'
i_string='define'
i_sunmath='undef'
i_sysaccess='undef'
i_sysdir='undef'
i_sysfile='undef'
i_sysfilio='define'
i_sysin='undef'
i_sysioctl='undef'
i_syslog='undef'
i_sysmman='undef'
i_sysmode='undef'
i_sysmount='undef'
i_sysndir='undef'
i_sysparam='undef'
i_sysresrc='undef'
i_syssecrt='undef'
i_sysselct='undef'
i_syssockio='undef'
i_sysstat='define'
i_sysstatfs='undef'
i_sysstatvfs='undef'
i_systime='undef'
i_systimek='undef'
i_systimes='undef'
i_systypes='define'
i_sysuio='undef'
i_sysun='undef'
i_sysutsname='undef'
i_sysvfs='undef'
i_syswait='undef'
i_termio='undef'
i_termios='undef'
i_time='define'
i_unistd='undef'
i_ustat='undef'
i_utime='define'
i_values='undef'
i_varargs='undef'
i_varhdr='varargs.h'
i_vfork='undef'
ignore_versioned_solibs=''
inc_version_list=''
inc_version_list_init='0'
incpath='C:\strawberry\c\include'
inews=''
initialinstalllocation=''
installarchlib='C:\strawberry\perl\lib'
installbin='C:\strawberry\perl\bin'
installhtml1dir=''
installhtml3dir=''
installhtmldir=''
installhtmlhelpdir=''
installman1dir=''
installman3dir=''
installprefix='C:\strawberry\perl'
installprefixexp='C:\strawberry\perl'
installprivlib='C:\strawberry\perl\lib'
installscript='C:\strawberry\perl\bin'
installsitearch='C:\strawberry\perl\site\lib'
installsitebin='C:\strawberry\perl\bin'
installsitehtml1dir=''
installsitehtml3dir=''
installsitelib='C:\strawberry\perl\site\lib'
installsiteman1dir=''
installsiteman3dir=''
installsitescript=''
installstyle='lib'
installusrbinperl='undef'
installvendorarch=''
installvendorbin=''
installvendorhtml1dir=''
installvendorhtml3dir=''
installvendorlib=''
installvendorman1dir=''
installvendorman3dir=''
installvendorscript=''
intsize='4'
issymlink=''
ivdformat='"ld"'
ivsize='4'
ivtype='long'
known_extensions='B Compress/Raw/Zlib Cwd DB_File Data/Dumper Devel/DProf Devel/PPPort Devel/Peek Digest/MD5 Digest/SHA Encode Fcntl File/Glob Filter/Util/Call GDBM_File Hash/Util Hash/Util/FieldHash I18N/Langinfo IO IPC/SysV List/Util MIME/Base64 Math/BigInt/FastCalc NDBM_File ODBM_File Opcode POSIX PerlIO/encoding PerlIO/scalar PerlIO/via SDBM_File Socket Storable Sys/Hostname Sys/Syslog Text/Soundex Time/HiRes Time/Piece Unicode/Normalize Win32 Win32API/File Win32CORE XS/APItest XS/Typemap attrs re threads threads/shared'
ksh=''
ld='g++'
lddlflags='-mdll -s -L"C:\strawberry\perl\lib\CORE" -L"C:\strawberry\c\lib"'
ldflags='-s -L"C:\strawberry\perl\lib\CORE" -L"C:\strawberry\c\lib"'
ldflags_uselargefiles=''
ldlibpthname=''
less='less'
lib_ext='.a'
libc='-lmsvcrt'
libperl='libperl510.a'
libpth='C:\strawberry\c\lib'
libs=' -lmsvcrt -lmoldname -lkernel32 -luser32 -lgdi32 -lwinspool -lcomdlg32 -ladvapi32 -lshell32 -lole32 -loleaut32 -lnetapi32 -luuid -lws2_32 -lmpr -lwinmm -lversion -lodbc32 -lodbccp32'
libsdirs=''
libsfiles=''
libsfound=''
libspath=''
libswanted='net socket inet nsl nm ndbm gdbm dbm db malloc dl dld ld sun m c cposix posix ndir dir crypt ucb bsd BSD PW x'
libswanted_uselargefiles='net socket inet nsl nm ndbm gdbm dbm db malloc dl dld ld sun m c cposix posix ndir dir crypt ucb bsd BSD PW x'
line='line'
lint=''
lkflags=''
ln=''
lns='copy'
localtime_r_proto='0'
locincpth='/usr/local/include /opt/local/include /usr/gnu/include /opt/gnu/include /usr/GNU/include /opt/GNU/include'
loclibpth='/usr/local/lib /opt/local/lib /usr/gnu/lib /opt/gnu/lib /usr/GNU/lib /opt/GNU/lib'
longdblsize='12'
longlongsize='8'
longsize='4'
lp=''
lpr=''
ls='dir'
lseeksize='8'
lseektype='long long'
mad='undef'
madlyh=''
madlyobj=''
madlysrc=''
mail=''
mailx=''
make='dmake'
make_set_make='#'
mallocobj='malloc.o'
mallocsrc='malloc.c'
malloctype='void *'
man1dir=''
man1direxp=''
man1ext='1'
man3dir=''
man3direxp=''
man3ext='3'
mips_type=''
mistrustnm=''
mkdir='mkdir'
mmaptype='void *'
modetype='mode_t'
more='more /e'
multiarch='undef'
mv=''
myarchname='MSWin32'
mydomain=''
myhostname=''
myuname=''
n='-n'
need_va_copy='undef'
netdb_hlen_type='int'
netdb_host_type='char *'
netdb_name_type='char *'
netdb_net_type='long'
nm='nm'
nm_opt=''
nm_so_opt=''
nonxs_ext='Compress/Zlib Errno IO_Compress_Base IO_Compress_Zlib'
nroff=''
nvEUformat='"E"'
nvFUformat='"F"'
nvGUformat='"G"'
nv_preserves_uv_bits='32'
nveformat='"e"'
nvfformat='"f"'
nvgformat='"g"'
nvsize='8'
nvtype='double'
o_nonblock='O_NONBLOCK'
obj_ext='.o'
old_pthread_create_joinable=''
optimize='-s -O2'
orderlib='false'
osname='MSWin32'
osvers='5.1'
otherlibdirs=''
package='perl5'
pager='more /e'
passcat=''
patchlevel='10'
path_sep=';'
perl5=''
perl='perl'
perl_patchlevel=''
perladmin=''
perllibs=' -lmsvcrt -lmoldname -lkernel32 -luser32 -lgdi32 -lwinspool -lcomdlg32 -ladvapi32 -lshell32 -lole32 -loleaut32 -lnetapi32 -luuid -lws2_32 -lmpr -lwinmm -lversion -lodbc32 -lodbccp32'
perlpath='C:\strawberry\perl\bin\perl.exe'
pg=''
phostname='hostname'
pidtype='int'
plibpth=''
pmake=''
pr=''
prefix='C:\strawberry\perl'
prefixexp='C:\strawberry\perl'
privlib='C:\strawberry\perl\lib'
privlibexp='C:\strawberry\perl\lib'
procselfexe=''
prototype='define'
ptrsize='4'
quadkind='5'
quadtype='long long'
randbits='15'
randfunc='rand'
random_r_proto='0'
randseedtype='unsigned'
ranlib='rem'
rd_nodata='-1'
readdir64_r_proto='0'
readdir_r_proto='0'
revision='5'
rm='del'
rmail=''
run=''
runnm='true'
sPRIEUldbl='"E"'
sPRIFUldbl='"F"'
sPRIGUldbl='"G"'
sPRIXU64='"lX"'
sPRId64='"ld"'
sPRIeldbl='"e"'
sPRIfldbl='"f"'
sPRIgldbl='"g"'
sPRIi64='"li"'
sPRIo64='"lo"'
sPRIu64='"lu"'
sPRIx64='"lx"'
sSCNfldbl='"f"'
sched_yield=''
scriptdir='C:\strawberry\perl\bin'
scriptdirexp='C:\strawberry\perl\bin'
sed='sed'
seedfunc='srand'
selectminbits='32'
selecttype='Perl_fd_set *'
sendmail='blat'
setgrent_r_proto='0'
sethostent_r_proto='0'
setlocale_r_proto='0'
setnetent_r_proto='0'
setprotoent_r_proto='0'
setpwent_r_proto='0'
setservent_r_proto='0'
sh='cmd /x /c'
shar=''
sharpbang='#!'
shmattype='void *'
shortsize='2'
shrpenv=''
shsharp='true'
sig_count='26'
sig_name='ZERO HUP INT QUIT ILL NUM05 NUM06 NUM07 FPE KILL NUM10 SEGV NUM12 PIPE ALRM TERM NUM16 NUM17 NUM18 NUM19 CHLD BREAK ABRT STOP NUM24 CONT CLD'
sig_name_init='"ZERO", "HUP", "INT", "QUIT", "ILL", "NUM05", "NUM06", "NUM07", "FPE", "KILL", "NUM10", "SEGV", "NUM12", "PIPE", "ALRM", "TERM", "NUM16", "NUM17", "NUM18", "NUM19", "CHLD", "BREAK", "ABRT", "STOP", "NUM24", "CONT", "CLD", 0'
sig_num='0 1 2 21 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 20'
sig_num_init='0, 1, 2, 21, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 20, 0'
sig_size='27'
signal_t='void'
sitearch='C:\strawberry\perl\site\lib'
sitearchexp='C:\strawberry\perl\site\lib'
sitebin='C:\strawberry\perl\site\bin'
sitebinexp='C:\strawberry\perl\site\bin'
sitehtml1dir=''
sitehtml1direxp=''
sitehtml3dir=''
sitehtml3direxp=''
sitelib='C:\strawberry\perl\site\lib'
sitelib_stem=''
sitelibexp='C:\strawberry\perl\site\lib'
siteman1dir=''
siteman1direxp=''
siteman3dir=''
siteman3direxp=''
siteprefix='C:\strawberry\perl\site'
siteprefixexp='C:\strawberry\perl\site'
sitescript=''
sitescriptexp=''
sizesize='4'
sizetype='size_t'
sleep=''
smail=''
so='dll'
sockethdr=''
socketlib=''
socksizetype='int'
sort='sort'
spackage='Perl5'
spitshell=''
srand48_r_proto='0'
srandom_r_proto='0'
src=''
ssizetype='int'
startperl='#!perl'
startsh='#!/bin/sh'
static_ext='Win32CORE'
stdchar='char'
stdio_base='((fp)->_base)'
stdio_bufsiz='((fp)->_cnt + (fp)->_ptr - (fp)->_base)'
stdio_cnt='((fp)->_cnt)'
stdio_filbuf=''
stdio_ptr='((fp)->_ptr)'
stdio_stream_array=''
strerror_r_proto='0'
strings='/usr/include/string.h'
submit=''
subversion=''
sysman=''
tail=''
tar=''
targetarch=''
tbl=''
tee=''
test=''
timeincl='/usr/include/sys/time.h '
timetype='time_t'
tmpnam_r_proto='0'
to=':'
touch='touch'
tr=''
trnl='\012'
troff=''
ttyname_r_proto='0'
u16size='2'
u16type='unsigned short'
u32size='4'
u32type='unsigned long'
u64size='8'
u64type='unsigned long long'
u8size='1'
u8type='unsigned char'
uidformat='"ld"'
uidsign='-1'
uidsize='4'
uidtype='uid_t'
uname='uname'
uniq='uniq'
uquadtype='unsigned long long'
use5005threads='undef'
use64bitall='undef'
use64bitint='undef'
usecrosscompile='undef'
usedl='define'
usefaststdio='undef'
useithreads='define'
uselargefiles='define'
uselongdouble='undef'
usemallocwrap='define'
usemorebits='undef'
usemultiplicity='define'
usemymalloc='n'
usenm='false'
useopcode='true'
useperlio='define'
useposix='true'
usereentrant='undef'
userelocatableinc='undef'
usesfio='false'
useshrplib='true'
usesitecustomize='undef'
usesocks='undef'
usethreads='define'
usevendorprefix='undef'
usevfork='false'
usrinc='/usr/include'
uuname=''
uvXUformat='"lX"'
uvoformat='"lo"'
uvsize='4'
uvtype='unsigned long'
uvuformat='"lu"'
uvxformat='"lx"'
vendorarch=''
vendorarchexp=''
vendorbin=''
vendorbinexp=''
vendorhtml1dir=''
vendorhtml1direxp=''
vendorhtml3dir=''
vendorhtml3direxp=''
vendorlib=''
vendorlib_stem=''
vendorlibexp=''
vendorman1dir=''
vendorman1direxp=''
vendorman3dir=''
vendorman3direxp=''
vendorprefix=''
vendorprefixexp=''
vendorscript=''
vendorscriptexp=''
version='5.10.0'
version_patchlevel_string='version 10 subversion 0'
versiononly='undef'
vi=''
voidflags='15'
xlibpth='/usr/lib/386 /lib/386'
yacc='yacc'
yaccflags=''
zcat=''
zip='zip'
!END!

my $i = 0;
foreach my $c (4,3,2) { $i |= ord($c); $i <<= 8 }
$i |= ord(1);
our $byteorder = join('', unpack('aaaa', pack('L!', $i)));
s/(byteorder=)(['"]).*?\2/$1$2$Config::byteorder$2/m;

my $config_sh_len = length $_;

our $Config_SH_expanded = "\n$_" . << 'EOVIRTUAL';
ccflags_nolargefiles=' -s -O2 -DWIN32 -DHAVE_DES_FCRYPT  -DPERL_IMPLICIT_CONTEXT -DPERL_IMPLICIT_SYS -fno-strict-aliasing -DPERL_MSVCRT_READFIX'
ldflags_nolargefiles='-s -L"C:\strawberry\perl\lib\CORE" -L"C:\strawberry\c\lib"'
libs_nolargefiles='-lmsvcrt -lmoldname -lkernel32 -luser32 -lgdi32 -lwinspool -lcomdlg32 -ladvapi32 -lshell32 -lole32 -loleaut32 -lnetapi32 -luuid -lws2_32 -lmpr -lwinmm -lversion -lodbc32 -lodbccp32'
libswanted_nolargefiles=''
EOVIRTUAL

# Search for it in the big string
sub fetch_string {
    my($self, $key) = @_;

    # We only have ' delimted.
    my $start = index($Config_SH_expanded, "\n$key=\'");
    # Start can never be -1 now, as we've rigged the long string we're
    # searching with an initial dummy newline.
    return undef if $start == -1;

    $start += length($key) + 3;

    my $value = substr($Config_SH_expanded, $start,
                       index($Config_SH_expanded, "'\n", $start)
		       - $start);
    # So we can say "if $Config{'foo'}".
    $value = undef if $value eq 'undef';
    $self->{$key} = $value; # cache it
}

my $prevpos = 0;

sub FIRSTKEY {
    $prevpos = 0;
    substr($Config_SH_expanded, 1, index($Config_SH_expanded, '=') - 1 );
}

sub NEXTKEY {
    my $pos = index($Config_SH_expanded, qq('\n), $prevpos) + 2;
    my $len = index($Config_SH_expanded, "=", $pos) - $pos;
    $prevpos = $pos;
    $len > 0 ? substr($Config_SH_expanded, $pos, $len) : undef;
}

sub EXISTS {
    return 1 if exists($_[0]->{$_[1]});

    return(index($Config_SH_expanded, "\n$_[1]='") != -1
          );
}

sub STORE  { die "\%Config::Config is read-only\n" }
*DELETE = \&STORE;
*CLEAR  = \&STORE;


sub config_sh {
    substr $Config_SH_expanded, 1, $config_sh_len;
}

sub config_re {
    my $re = shift;
    return map { chomp; $_ } grep eval{ /^(?:$re)=/ }, split /^/,
    $Config_SH_expanded;
}

sub config_vars {
    # implements -V:cfgvar option (see perlrun -V:)
    foreach (@_) {
	# find optional leading, trailing colons; and query-spec
	my ($notag,$qry,$lncont) = m/^(:)?(.*?)(:)?$/;	# flags fore and aft, 
	# map colon-flags to print decorations
	my $prfx = $notag ? '': "$qry=";		# tag-prefix for print
	my $lnend = $lncont ? ' ' : ";\n";		# line ending for print

	# all config-vars are by definition \w only, any \W means regex
	if ($qry =~ /\W/) {
	    my @matches = config_re($qry);
	    print map "$_$lnend", @matches ? @matches : "$qry: not found"		if !$notag;
	    print map { s/\w+=//; "$_$lnend" } @matches ? @matches : "$qry: not found"	if  $notag;
	} else {
	    my $v = (exists $Config::Config{$qry}) ? $Config::Config{$qry}
						   : 'UNKNOWN';
	    $v = 'undef' unless defined $v;
	    print "${prfx}'${v}'$lnend";
	}
    }
}

# Called by the real AUTOLOAD
sub launcher {
    undef &AUTOLOAD;
    goto \&$Config::AUTOLOAD;
}

1;
