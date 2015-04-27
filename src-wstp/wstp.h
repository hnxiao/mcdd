/*************************************************************************

        Copyright 1986 through 2014 by Wolfram Research Inc.
        All rights reserved

*************************************************************************/

#ifndef _WSTP_H
#define _WSTP_H

#ifndef _WSVERS_H
#define _WSVERS_H

#ifndef _WSPLATFM_H
#define _WSPLATFM_H

#if ! WINDOWS_WSTP && ! UNIX_WSTP
#	define	UNIX_WSTP	1
#endif

#if defined(__GNUC__)
#if defined(__GNUC_PATCHLEVEL__)
#define GCC_WSTP_VERSION (__GNUC__ * 10000 \
						+ __GNUC_MINOR__ * 100 \
						+ __GNUC_PATCHLEVEL__)
#else
#define GCC_WSTP_VERSION (__GNUC__ * 10000 \
						+ __GNUC_MINOR__ * 100)
#endif
#endif

#if WINDOWS_WSTP
#if defined(WIN64) || defined(__WIN64__) || defined(_WIN64)
#define WIN64_WSTP 1
#if( _M_IX86 || __i386 || __i386__ || i386)
#define I86_WIN32_WSTP 1
#endif
#elif defined(WIN32) || defined(__WIN32__) || defined(__NT__) || defined(_WIN32)
#define WIN32_WSTP 1
#if( _M_IX86 || __i386 || __i386__ || i386)
#define I86_WIN32_WSTP 1
#elif _M_ALPHA || __alpha || __alpha__ || alpha
#define ALPHA_WIN32_WSTP 1
#else
#endif
#elif defined(__CYGWIN__)
#define CYGWIN_WSTP 1
#else
#define WIN16_WSTP 1
#endif
#elif UNIX_WSTP
#if (__sun || __sun__ || sun) && !defined(SUN_WSTP)
#define SUN_WSTP 1
#if __sparcv9 || __sparcv9__
#if __SUNPRO_C >= 0x590
#define STUDIO12_U64_SOLARIS_WSTP 1
#else
#define U64_SOLARIS_WSTP 1
#endif
#elif __SVR4 || __svr4__
#define U32_SOLARIS_WSTP 1
#else
#define SUNOS_WSTP 1
#endif
#if __sparc || __sparc__ || sparc
#define SPARC_SUN_WSTP 1
#elif __i386 || __i386__ || i386
#define I386_SOLARIS_WSTP 1
#elif __x86_64 || __x86_64__ || x86_64
#define X86_64_SOLARIS_WSTP 1
#else
			unknown platform
#endif
#elif (__MACH || __MACH__ || MACH) && !defined(DARWIN_WSTP)
#define DARWIN_WSTP 1
#if __ppc || __ppc__ || ppc
#define POWERPC_DARWIN_WSTP 1
#define PPC_DARWIN_WSTP 1
#elif __ppc64 || __ppc64__ || ppc64
#define POWERPC_DARWIN_WSTP 1
#define PPC64_DARWIN_WSTP 1
#elif __i386 || __i386__ || i386
#define INTEL_DARWIN_WSTP 1
#define X86_DARWIN_WSTP 1
#elif __x86_64 || __x86_64__ || x86_64
#define INTEL_DARWIN_WSTP 1
#define X86_64_DARWIN_WSTP 1
#elif __arm__
#define ARM_DARWIN_WSTP 1
#elif __arm64__
#define ARM64_DARWIN_WSTP 1
#else
			not yet implemented
#endif

#if __DARWIN_UNIX03
#define DARWIN_WSTP_UNIX03 1
#endif /* __DARWIN_UNIX03 */

#if __clang__
#define CLANG_WSTP 1
#endif
#elif (__linux || __linux__ || linux) && !defined(LINUX_WSTP)
#define LINUX_WSTP 1
#if __x86_64 || __x86_64__ || x86_64
#define X86_64_LINUX_WSTP 1
#elif __ia64 || __ia64__ || ia64
#define IA64_LINUX_WSTP 1
#elif __i386 || __i386__ || i386
#define I86_LINUX_WSTP 1
#elif __PPC || __PPC__ || PPC
#define PPC_LINUX_WSTP 1
#elif __alpha || __alpha__ || alpha
#define AXP_LINUX_WSTP 1
#elif __ANDROID || __ANDROID__ || ANDROID
#define ANDROID_LINUX_WSTP 1
#if __arm || __arm__ || arm
#define ANDROID_ARM_LINUX_WSTP 1
#endif
#elif __arm || __arm__ || arm
#define ARM_LINUX_WSTP 1
#else
			not yet implemented
#endif
#elif (__osf || __osf__ || osf || OSF1) && !defined(DIGITAL_WSTP)
#define DIGITAL_WSTP 1
#if __alpha || __alpha__ || alpha
#define ALPHA_DIGITAL_WSTP 1
#else
			unknown platform
#endif
#elif (_AIX || _IBMR2 || __xlC__) && !defined(AIX_WSTP)
#define AIX_WSTP 1
#if __64BIT__
#define A64_AIX_WSTP 1
#endif
#elif (__sgi || __sgi__ || sgi || mips) && !defined(IRIX_WSTP)
#define IRIX_WSTP 1
#if _MIPS_SZLONG == 32
#define N32_IRIX_WSTP 1
#elif _MIPS_SZLONG == 64
#define M64_IRIX_WSTP 1
#else
			not yet implemented
#endif
#elif (hpux || __hpux) && !defined(HPUX_WSTP)
#define HPUX_WSTP 1
#if __LP64__
#define LP64_HPUX_WSTP 1
#endif
#elif (M_I386 || _SCO_DS || SCO) && !defined(SCO_WSTP)
#define SCO_WSTP 1
#elif (__NetBSD__) && !defined(NETBSD_WSTP)
#define NETBSD_WSTP 1
#elif (__FreeBSD__) && !defined(FREEBSD_WSTP)
#define FREEBSD_WSTP 1
#elif (bsdi || __bsdi__) && !defined(BSDI_WSTP)
#define BSDI_WSTP 1
#else
#endif
#endif



#ifndef NO_GLOBAL_DATA
#define NO_GLOBAL_DATA 0
#endif

#if WINDOWS_WSTP || __i386 || __i386__ || i386 || _M_IX86 || __x86_64 || __x86_64__ || x86_64 || __ia64 || __ia64__ || ia64 ||  __alpha || __alpha__ || alpha || __arm__
#define LITTLEENDIAN_NUMERIC_TYPES 1
#else
#define BIGENDIAN_NUMERIC_TYPES 1
#endif

#endif /* _WSPLATFM_H */



#if WSINTERFACE <= 3
#ifndef WSVERSION
#define WSVERSION 4
#endif
#endif

#if WSINTERFACE >= 4
#ifndef WSVERSION
#define WSVERSION 5
#endif
#endif

#if !OLD_VERSIONING


/*
 * WSTP adopts a simple versioning strategy that can be adapted to many
 * compile-time and run-time environments.  In particular, it is amenable to
 * the various shared library facilities in use.  (Although certain of these
 * facilities provide more sophisticated mechanisms than are required by the
 * following simple strategy.)
 * 
 * WSTP evolves by improving its implementation and by improving its
 * interface.  The values of WSREVISION or WSINTERFACE defined here are 
 * incremented whenever an improvement is made and released.
 * 
 * WSREVISION is the current revision number. It is incremented every time 
 * a change is made to the source and WSTP is rebuilt and distributed
 * on any platform.  (Bug fixes, optimizations, or other improvements
 * transparent to the interface increment only this number.)
 * 
 * WSINTERFACE is a name for a documented interface to WSTP.  This
 * number is incremented whenever a named constant or function is added,
 * removed, or its behavior is changed in a way that could break existing
 * correct* client programs.  It is expected that the interface to WSTP
 * is improved over time so that implemenations of higher numbered 
 * interfaces are more complete or more convenient to use for writing
 * effective client programs.  In particular, a specific interface provides
 * all the useful functionality of an earlier interface.
 * 
 *     *(It is possible that an incorrect WSTP program still works
 *     because it relies on some undocumented detail of a particular
 *     revision.  It may not always be possible to change the interface
 *     number when such a detail changes.  For example, one program may
 *     be relying on a bug in WSTP that a great many other programs
 *     need fixed.  In this case, we would likely choose to potentially
 *     break the incorrect program in order to fix the correct programs
 *     by incrementing the revision number leaving the interface number
 *     unchanged.  It is possible to bind to a particular revision of a
 *     WSTP interface if that is important for some programs.  One
 *     could use a statically linked version of the library, make use of
 *     the search algorithm used by the runtime loader, or dynamically
 *     load the WSTP library manually.)
 * 
 * 
 * If a distributed WSTP implmentation were labeled with its revision 
 * and interface numbers in dotted notation so that, say, WS.1.6 means the
 * sixth revision of interface one, then the following may represent the
 * distribution history of WSTP.
 * 
 *     first distribution
 *         WS.1.5   (Perhaps earlier revisions were never
 *                   distributed for this platform.)
 * 
 *     second distribution
 *         WS.1.6   (Bug fixes or other improvements were
 *                   made that don't affect the interface.)
 * 
 *     third distribution
 *         WS.2.7   (Perhaps some new functions were added.)
 *         
 *         WS.1.7   (And improvements were made that don't
 *                   affect the old interface.)
 * 
 *     fourth distribution
 *         WS.3.8   (Perhaps the return values of an existing
 *                   function changed.)
 *         WS.2.8   (Revision 8 also adds improvements transparent
 *                   to interface 2.)
 *         WS.1.7   (Clients of interface 1 see no improvements
 *                   in this eighth revision.)
 * 
 * Note that the distribution history may not be the same on different
 * platforms.  But revision numbers represent a named body of source code
 * across all platforms.
 * 
 * The mechanism for deploying this strategy differs between platforms
 * because of differing platform-specific facilities and conventions.
 * The interface and revision numbers may form part of the filename of
 * the WSTP library, or they may not.  This information is always
 * available in some conventional form so that it is easy and natural for
 * client programs to bind with and use the best available implementation
 * of a particular WSTP interface.  The details are described in the
 * WSTP Developer's Guide for each platform.
 */

#define WSREVISION 25
#define WSMATHVERSION 10.0.0

#ifndef WSCREATIONID
#define WSCREATIONID 114411
#endif

#define WSAPI1REVISION 1   /* the first revision to support interface 1 */
#define WSAPI2REVISION 6   /* the first revision to support interface 2 */
#define WSAPI3REVISION 16  /* the first revision to support interface 3 */
#define WSAPI4REVISION 25  /* the first revision to support interface 4 */


#ifndef WSINTERFACE
#define WSINTERFACE 4
#define WSAPIREVISION WSAPI4REVISION

		/*
		 * Interface 4 adds the following exported functions:
         *
         * WSCreateLinkWithExternalProtocol
         * WSDoNotHandleSignalParameter
         * WSEnableLinkLock
         * WSFilterArgv
         * WSGetAvailableLinkProtocolNames
         * WSGetInteger8
         * WSGetInteger8Array
         * WSGetInteger8ArrayData
         * WSGetInteger8List
         * WSGetLinksFromEnvironment
         * WSGetNumberAsByteString
         * WSGetNumberAsString
         * WSGetNumberAsUCS2String
         * WSGetNumberAsUTF16String
         * WSGetNumberAsUTF32String
         * WSGetNumberAsUTF8String
         * WSHandleSignal
         * WSIsLinkLoopback
         * WSLinkEnvironment
         * WSLogFileNameForLink
         * WSLogStreamToFile
         * WSLowLevelDeviceName
         * WSPutInteger8
         * WSPutInteger8Array
         * WSPutInteger8ArrayData
         * WSPutInteger8List
         * WSPutRealNumberAsUCS2String
         * WSPutRealNumberAsUTF16String
         * WSPutRealNumberAsUTF32String
         * WSPutRealNumberAsUTF8String
         * WSReleaseInteger8Array
         * WSReleaseInteger8List
         * WSReleaseLinkProtocolNames
         * WSReleaseLinksFromEnvironment
         * WSReleaseLogFileNameForLink
         * WSReleaseLowLevelDeviceName
         * WSReleaseParameters
         * WSSetThreadSafeLinksParameter
         * WSStopHandlingSignal
         * WSStopLoggingStream
         * WSValid
         * WSWaitForLinkActivity
         * WSWaitForLinkActivityWithCallback
         *
         * Interface 4 removes the following API functions.
         *
         * WSCreate0
         * WSDestroy
         * WSDoNotHandleSignalParameter0
         * WSFilterArgv0
         * WSGetByteString0
         * WSGetLinkedEnvIDString0
         * WSGetString0
         * WSGetUCS2String0
         * WSGetUTF16String0
         * WSGetUTF32String0
         * WSGetUTF8String0
         * WSGetUnicodeString0
         * WSHandleSignal0
         * WSLoopbackOpen0
         * WSMake
         * WSPutRealByteString0
         * WSPutRealUCS2String0
         * WSPutRealUTF16String0
         * WSPutRealUTF32String0
         * WSPutRealUTF8String0
         * WSPutRealUnicodeString0
         * WSSetEnvIDString0
         * WSSetMessageHandler0
         * WSSetSignalHandler0
         * WSSetYieldFunction0
         * WSStopHandlingSignal0
         * WSUnsetSignalHandler0
         * WSValid0
         * WSVersionNumber0
         *
		 */


                /*
                 * Interface 3 adds the following exported functions:
                 *      WSClearAllSymbolReplacements
                 *      WSClearSymbolReplacement
                 *      WSConvertUCS2String
                 *      WSConvertUCS2StringNL
                 *      WSConvertUTF8String
                 *      WSConvertUTF8StringNL
                 *      WSConvertUTF16String
                 *      WSConvertUTF16StringNL
                 *      WSConvertUTF32String
                 *      WSConvertUTF32StringNL
                 *      WSEnvironmentData
                 *      WSGetDomainNameList
                 *      WSGetInteger16
                 *      WSGetInteger16Array
                 *      WSGetInteger16ArrayData
                 *      WSGetInteger16List
                 *      WSGetInteger32
                 *      WSGetInteger32Array
                 *      WSGetInteger32ArrayData
                 *      WSGetInteger32List
                 *      WSGetInteger64
                 *      WSGetInteger64Array
                 *      WSGetInteger64ArrayData
                 *      WSGetInteger64List
                 *      WSGetLinkedEnvIDString
                 *      WSGetMessageHandler
                 *      WSGetNetworkAddressList
                 *      WSGetReal128
                 *      WSGetReal128Array
                 *      WSGetReal128ArrayData
                 *      WSGetReal128List
                 *      WSGetReal32
                 *      WSGetReal32Array
                 *      WSGetReal32ArrayData
                 *      WSGetReal32List
                 *      WSGetReal64
                 *      WSGetReal64Array
                 *      WSGetReal64ArrayData
                 *      WSGetReal64List
                 *      WSGetUCS2Characters
                 *      WSGetUCS2String
                 *      WSGetUCS2Symbol
                 *      WSGetUTF16Characters
                 *      WSGetUTF16String
                 *      WSGetUTF16Symbol
                 *      WSGetUTF32Characters
                 *      WSGetUTF32String
                 *      WSGetUTF32Symbol
                 *      WSGetUTF8Characters
                 *      WSGetUTF8String
                 *      WSGetUTF8Symbol
                 *      WSGetYieldFunction
                 *      WSLinkName
                 *      WSOldConvertUCS2String
                 *      WSOpenArgcArgv
                 *      WSPutInteger16
                 *      WSPutInteger16Array
                 *      WSPutInteger16ArrayData
                 *      WSPutInteger16List
                 *      WSPutInteger32
                 *      WSPutInteger32Array
                 *      WSPutInteger32ArrayData
                 *      WSPutInteger32List
                 *      WSPutInteger64
                 *      WSPutInteger64Array
                 *      WSPutInteger64ArrayData
                 *      WSPutInteger64List
                 *      WSPutMessageWithArg
                 *      WSPutReal128
                 *      WSPutReal128Array
                 *      WSPutReal128ArrayData
                 *      WSPutReal128List
                 *      WSPutReal32
                 *      WSPutReal32Array
                 *      WSPutReal32ArrayData
                 *      WSPutReal32List
                 *      WSPutReal64
                 *      WSPutReal64Array
                 *      WSPutReal64ArrayData
                 *      WSPutReal64List
                 *      WSPutUCS2Characters
                 *      WSPutUCS2String
                 *      WSPutUCS2Symbol
                 *      WSPutUTF16Characters
                 *      WSPutUTF16String
                 *      WSPutUTF16Symbol
                 *      WSPutUTF32Characters
                 *      WSPutUTF32String
                 *      WSPutUTF32Symbol
                 *      WSPutUTF8Characters
                 *      WSPutUTF8String
                 *      WSPutUTF8Symbol
                 *      WSReadyParallel
                 *      WSReleaseBinaryNumberArray
                 *      WSReleaseByteArray
                 *      WSReleaseByteString
                 *      WSReleaseByteSymbol
                 *      WSReleaseDomainNameList
                 *      WSReleaseInteger16Array
                 *      WSReleaseInteger16List
                 *      WSReleaseInteger32Array
                 *      WSReleaseInteger32List
                 *      WSReleaseInteger64Array
                 *      WSReleaseInteger64List
                 *      WSReleaseNetworkAddressList
                 *      WSReleaseReal128Array
                 *      WSReleaseReal128List
                 *      WSReleaseReal32Array
                 *      WSReleaseReal32List
                 *      WSReleaseReal64Array
                 *      WSReleaseReal64List
                 *      WSReleaseString
                 *      WSReleaseSymbol
                 *      WSReleaseUCS2String
                 *      WSReleaseUCS2Symbol
                 *      WSReleaseUTF16String
                 *      WSReleaseUTF16Symbol
                 *      WSReleaseUTF32String
                 *      WSReleaseUTF32Symbol
                 *      WSReleaseUTF8String
                 *      WSReleaseUTF8Symbol
                 *      WSSetEncodingParameter
                 *      WSSetEnvIDString
                 *      WSSetEnvironmentData
                 *      WSSetSignalHandler
                 *      WSSetSignalHandlerFromFunction
                 *      WSSetSymbolReplacement
                 *      WSTestHead
                 *      WSUnsetSignalHandler
                 *      WSVersionNumbers
                 *      
                 *      Interface 3 removes all the special WSTP types such as kcharp_ct, long_st, etc.
                 *      The API functions now all take standard C data types.
                 */

       		/*
                 * Interface 2 adds the following exported functions:
                 *      WSGetBinaryNumberArray0
                 *      WSTransfer0
                 *      WSNextCharacter0
                 * And, for WINDOWS_WSTP, some numerics constants
                 * were changed in a way that causes WSGetRawType to return
                 * different values.
                 *
                 */
#else
#if WSINTERFACE == 1
#define WSAPIREVISION WSAPI1REVISION
#elif WSINTERFACE == 2
#define WSAPIREVISION WSAPI2REVISION
#elif WSINTERFACE == 3
#define WSAPIREVISION WSAPI3REVISION
#elif WSINTERFACE == 4
#define WSAPIREVISION WSAPI4REVISION
#else
/*\ syntax\ error\ */\ )
#endif
#endif


/* It may be possible for an implementation of one WSTP interface to
 * fully support an earlier interface.  WSNewParameters() may succeed when
 * passed an interface number less than the value of WSAPIREVISION when the
 * library was built.  This would happen, if the newer interface is a proper
 * superset of the older interface, or if the implementation can adjust its
 * behavior at runtime to conform to the older requested interface.
 */

#ifndef WSOLDDEFINITION
#if WINDOWS_WSTP
#if WSINTERFACE == 1
#define WSOLDDEFINITION WSAPI1REVISION
#elif WSINTERFACE == 2
#define WSOLDDEFINITION WSAPI2REVISION
#elif WSINTERFACE == 3
#define WSOLDDEFINITION WSAPI2REVISION
#elif WSINTERFACE == 4
#define WSOLDDEFINITION WSAPI3REVISION
#else
/*\ syntax\ error\ */\ )
#endif
#else
#define WSOLDDEFINITION WSAPI1REVISION
#endif
#endif




#else
/*\ syntax\ error\ */\ )
#endif

#endif /* _WSVERS_H */




#ifndef WS_EXTERN_C

#if defined(__cplusplus)
#define WS_C "C"
#define WS_EXTERN_C extern "C" {
#define WS_END_EXTERN_C }
#else
#define WS_C
#define WS_EXTERN_C
#define WS_END_EXTERN_C
#endif

#endif






#if WINDOWS_WSTP && (MPREP_REVISION || !defined(APIENTRY) || !defined(FAR))

#if defined(WIN32_LEAN_AND_MEAN) && defined(WIN32_EXTRA_LEAN)
#include <windows.h>
#elif defined( WIN32_LEAN_AND_MEAN)
#define WIN32_EXTRA_LEAN
#include <windows.h>
#undef WIN32_EXTRA_LEAN
#elif defined( WIN32_EXTRA_LEAN)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#else
#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_EXTRA_LEAN
#undef WIN32_LEAN_AND_MEAN
#endif



#endif


#ifndef _WSDEVICE_H
#define _WSDEVICE_H



#ifndef P


#ifndef WSPROTOTYPES
#define WSPROTOTYPES 1
#endif

#if WSPROTOTYPES || __STDC__ || defined(__cplusplus) || ! UNIX_WSTP
#define P(s) s
#undef WSPROTOTYPES
#define WSPROTOTYPES 1
#else
#define P(s) ()
#undef WSPROTOTYPES
#define WSPROTOTYPES 0
#endif
#endif



#ifndef _WSFAR_H
#define _WSFAR_H

#ifndef FAR


#if WINDOWS_WSTP

#ifndef FAR
/*\ syntax\ error\ */\ )
#endif
#else
#define FAR
#endif


#endif


#ifndef WSHUGE
#if WINDOWS_WSTP && !(WIN32_WSTP || WIN64_WSTP)
#define WSHUGE huge
#else
#define WSHUGE
#endif
#endif

#endif /* _WSFAR_H */



#ifndef _WSTYPES_H
#define _WSTYPES_H






#ifndef _WSBASICTYPES_H
#define _WSBASICTYPES_H




#ifndef _WSINT64_H
#define _WSINT64_H


#if WIN64_WSTP || X86_64_SOLARIS_WSTP || IA64_LINUX_WSTP || X86_64_LINUX_WSTP || A64_AIX_WSTP || M64_IRIX_WSTP || LP64_HPUX_WSTP || PPC64_DARWIN_WSTP || X86_64_DARWIN_WSTP || DIGITAL_WSTP || U64_SOLARIS_WSTP || STUDIO12_U64_SOLARIS_WSTP
#define WS64BIT_WSTP 1
#endif

#endif /* WSINT64_H */




#if WS64BIT_WSTP
#define WS_SMALLEST_SIGNED_64BIT        -9223372036854775807L - 1
#define WS_LARGEST_SIGNED_64BIT         9223372036854775807L
#define WS_LARGEST_UNSIGNED_64BIT       18446744073709551615UL
#else
#define WS_SMALLEST_SIGNED_64BIT        -9223372036854775807LL - 1
#define WS_LARGEST_SIGNED_64BIT         9223372036854775807LL
#define WS_LARGEST_UNSIGNED_64BIT       18446744073709551615ULL
#endif /* WS64BIT_WSTP */

#ifndef NO_INT64_STRUCT
#define NO_INT64_STRUCT

#if WS64BIT_WSTP
typedef unsigned int wint;
#else
typedef unsigned long wint;
#endif /* WS64BIT_WSTP */

#if LITTLEENDIAN_NUMERIC_TYPES
typedef struct _wint{
	wint low, hi;
} wint64;
#else
typedef struct _wint{
	wint hi, low;
} wint64;
#endif /* LITTLEENDIAN_NUMERIC_TYPES */


#endif /* NO_INT64_STRUCT */


#if WSINTERFACE >= 4
#if WIN64_WSTP || WIN32_WSTP

typedef long wslong32;
typedef unsigned long wsulong32;


typedef __int64 wsint64;
typedef unsigned __int64 wsuint64;

#if WIN64_WSTP
typedef wsint64 wsbigint;
typedef wsuint64 wsbiguint;
#else
typedef wslong32 wsbigint;
typedef wsulong32 wsbiguint;
#endif


#elif X86_64_SOLARIS_WSTP || IA64_LINUX_WSTP || X86_64_LINUX_WSTP || A64_AIX_WSTP || M64_IRIX_WSTP || LP64_HPUX_WSTP || PPC64_DARWIN_WSTP || X86_64_DARWIN_WSTP || DIGITAL_WSTP || U64_SOLARIS_WSTP || STUDIO12_U64_SOLARIS_WSTP

typedef int wslong32;
typedef unsigned int wsulong32;


typedef long wsint64;
typedef unsigned long wsuint64;

typedef wsint64 wsbigint;
typedef wsuint64 wsbiguint;

#else /* All other 32 bit platforms */

typedef long wslong32;
typedef unsigned long wsulong32;


typedef long long wsint64;
typedef unsigned long long wsuint64;

typedef wslong32 wsbigint;
typedef wsulong32 wsbiguint;

#endif /* WIN64_WSTP || WIN32_WSTP */



#endif /* WSINTERFACE >= 4 */




#if WSINTERFACE < 4

#if WIN64_WSTP


typedef __int64 wsint64;
typedef unsigned __int64 wsuint64;

#elif X86_64_SOLARIS_WSTP || IA64_LINUX_WSTP || X86_64_LINUX_WSTP || A64_AIX_WSTP || M64_IRIX_WSTP || LP64_HPUX_WSTP || PPC64_DARWIN_WSTP || X86_64_DARWIN_WSTP || DIGITAL_WSTP || U64_SOLARIS_WSTP || STUDIO12_U64_SOLARIS_WSTP


typedef long wsint64;
typedef unsigned long wsuint64;

#else /* All other 32 bit platforms */


#ifndef NO_WSINT64_STRUCT
#define NO_WSINT64_STRUCT

typedef wint64 wsint64;
typedef wint64 wsuint64;

#endif /* NO_WSINT64_STRUCT */

#endif


#endif /* WSINTERFACE < 4 */


#endif /* _WSBASICTYPES_H */




#if WINDOWS_WSTP

#ifndef	APIENTRY
#define APIENTRY far pascal
#endif
#ifndef CALLBACK
#define CALLBACK APIENTRY
#endif
#if (WIN32_WSTP || WIN64_WSTP || CYGWIN_WSTP)
#define WSEXPORT
#else
#define WSEXPORT __export
#endif
#define WSCB APIENTRY WSEXPORT
#define WSAPI APIENTRY

#else
#define WSCB
#define WSAPI
#define WSEXPORT
#endif

#if defined(LINUX_MATHLINK)
#define WSATTR __attribute__ ((visibility("default")))
#else
#define WSATTR
#endif

#define WSAPI_ WSAPI

#ifndef WSDEFN
#define WSDEFN( rtype, name, params) extern WSATTR rtype WSAPI WSEXPORT name params
#endif

#ifndef WSDECL
#define WSDECL( rtype, name, params) extern rtype WSAPI name P(params)
#endif

#ifndef WS_DEFN
#define WS_DEFN( rtype, name, params) extern rtype WSAPI_ WSEXPORT name params
#endif
#ifndef WS_DECL
#define WS_DECL( rtype, name, params) extern WS_C rtype WSAPI_ name P(params)
#endif



#ifndef WSCBPROC
#define WSCBPROC( rtype, name, params) typedef rtype (WSCB * name) P(params)
#endif
#ifndef WSCBDECL
#define WSCBDECL( rtype, name, params) extern rtype WSCB name P(params)
#endif
#ifndef WSCBDEFN
#define WSCBDEFN( rtype, name, params) extern rtype WSCB name params
#endif




#ifndef WSDPROC
#define WSDPROC WSCBPROC
#endif
#ifndef WSDDECL
#define WSDDECL WSCBDECL
#endif
#ifndef WSDDEFN
#define WSDDEFN WSCBDEFN
#endif




#ifndef WSTPROC
#define WSTPROC WSCBPROC
#endif
#ifndef WSTDECL
#define WSTDECL WSCBDECL
#endif
#ifndef WSTDEFN
#define WSTDEFN WSCBDEFN
#endif


#ifndef WSNPROC
#define WSNPROC WSCBPROC
#endif
#ifndef WSNDECL
#define WSNDECL WSCBDECL
#endif
#ifndef WSNDEFN
#define WSNDEFN WSCBDEFN
#endif


#ifndef WSAPROC
#define WSAPROC WSCBPROC
#endif
#ifndef WSADECL
#define WSADECL WSCBDECL
#endif
#ifndef WSADEFN
#define WSADEFN WSCBDEFN
#endif
#ifndef WSFPROC
#define WSFPROC WSCBPROC
#endif
#ifndef WSFDECL
#define WSFDECL WSCBDECL
#endif
#ifndef WSFDEFN
#define WSFDEFN WSCBDEFN
#endif
#ifndef WSRAPROC
#define WSRAPROC WSCBPROC
#endif
#ifndef WSRADECL
#define WSRADECL WSCBDECL
#endif
#ifndef WSRADEFN
#define WSRADEFN WSCBDEFN
#endif


#ifndef WSYPROC
#define WSYPROC WSCBPROC
#endif
#ifndef WSYDECL
#define WSYDECL WSCBDECL
#endif
#ifndef WSYDEFN
#define WSYDEFN WSCBDEFN
#endif
#ifndef WSMPROC
#define WSMPROC WSCBPROC
#endif
#ifndef WSMDECL
#define WSMDECL WSCBDECL
#endif
#ifndef WSMDEFN
#define WSMDEFN WSCBDEFN
#endif


#ifndef WSUPROC
#define WSUPROC WSCBPROC
#endif
#ifndef WSUDECL
#define WSUDECL WSCBDECL
#endif
#ifndef WSUDEFN
#define WSUDEFN WSCBDEFN
#endif


#ifndef WSBPROC
#define WSBPROC WSCBPROC
#endif
#ifndef WSBDECL
#define WSBDECL WSCBDECL
#endif
#ifndef WSBDEFN
#define WSBDEFN WSCBDEFN
#endif

#ifndef WSDMPROC
#define WSDMPROC WSCBPROC
#endif
#ifndef WSDMDECL
#define WSDMDECL WSCBDECL
#endif
#ifndef WSDMDEFN
#define WSDMDEFN WSCBDEFN
#endif

#if WSINTERFACE >= 4
#ifndef WSWPROC
#define WSWPROC WSCBPROC
#endif
#endif

#ifndef __uint_ct__
#define __uint_ct__ unsigned int
#endif
#ifndef __int_ct__
#define __int_ct__ int
#endif


typedef unsigned char        uchar_ct;
typedef uchar_ct       FAR * ucharp_ct;
typedef ucharp_ct      FAR * ucharpp_ct;
typedef ucharpp_ct     FAR * ucharppp_ct;
typedef unsigned short       ushort_ct;
typedef ushort_ct      FAR * ushortp_ct;
typedef ushortp_ct     FAR * ushortpp_ct;
typedef ushortpp_ct    FAR * ushortppp_ct;
typedef __uint_ct__          uint_ct;
typedef __uint_ct__    FAR * uintp_ct;
typedef uintp_ct       FAR * uintpp_ct;
typedef __int_ct__           int_ct;
typedef void           FAR * voidp_ct;
typedef voidp_ct       FAR * voidpp_ct;
typedef char           FAR * charp_ct;
typedef charp_ct       FAR * charpp_ct;
typedef charpp_ct      FAR * charppp_ct;
typedef long                 long_ct;
typedef long_ct        FAR * longp_ct;
typedef longp_ct       FAR * longpp_ct;
typedef unsigned long        ulong_ct;
typedef ulong_ct       FAR * ulongp_ct;




#ifndef WSCONST
#if WSPROTOTYPES
#define WSCONST const
#else
#define WSCONST
#endif
#endif

typedef WSCONST unsigned short FAR * kushortp_ct;
typedef WSCONST unsigned short FAR * FAR * kushortpp_ct;
typedef WSCONST unsigned int FAR * kuintp_ct;
typedef WSCONST unsigned int FAR * FAR * kuintpp_ct;
typedef WSCONST unsigned char FAR * kucharp_ct;
typedef WSCONST unsigned char FAR * FAR * kucharpp_ct;
typedef WSCONST char FAR * kcharp_ct;
typedef WSCONST char FAR * FAR * kcharpp_ct;
typedef WSCONST void FAR * kvoidp_ct;


typedef void FAR * WSPointer;

#ifndef __WSENVPARAM__
	typedef void FAR * WSENVPARAM;
	typedef WSENVPARAM WSEnvironmentParameter;
#define __WSENVPARAM__ 
#endif

#ifndef __WSENV__
	typedef struct ws_environment FAR *WSENV;
	typedef WSENV WSEnvironment;
#define __WSENV__
#endif

#ifndef __WSLINK__
	typedef struct WSLink FAR *WSLINK;
#define __WSLINK__
#endif

#ifndef __WSMARK__
	typedef struct WSLinkMark FAR *WSMARK;
	typedef WSMARK WSINKMark;
#define __WSMARK__
#endif

#ifndef __wsapi_token__
#define __wsapi_token__ int_ct
#endif
typedef __wsapi_token__   wsapi_token;


typedef unsigned long      wsapi__token;
typedef wsapi__token FAR * wsapi__tokenp;

#ifndef __wsapi_packet__
#define __wsapi_packet__ int_ct
#endif
typedef __wsapi_packet__  wsapi_packet;


typedef long wsapi_error;
typedef long wsapi__error;

typedef long  long_st;
typedef longp_ct longp_st;
typedef longp_st* longpp_st;

typedef long long_et;


#ifndef __wsapi_result__
#define __wsapi_result__ int_ct
#endif
typedef __wsapi_result__ wsapi_result;


#define WSSUCCESS (1) /*bugcheck:  this stuff doesnt belong where it can be seen at WSAPI_ layer */
#define WSFAILURE (0)

WS_EXTERN_C

#if WINDOWS_WSTP
typedef int (CALLBACK *__WSProcPtr__)();
#else
typedef long (* __WSProcPtr__)(void);
#endif

WS_END_EXTERN_C

#endif /* _WSTYPES_H */




#if WINDOWS_WSTP

#ifndef	APIENTRY
#define	APIENTRY far pascal
#endif
#define WSBN APIENTRY /* bottleneck function: upper layer calls lower layer */
#else
#define WSBN
#endif

#define BN WSBN



WS_EXTERN_C



typedef void FAR * dev_voidp;
typedef dev_voidp dev_type;
typedef dev_type FAR * dev_typep;
typedef long devproc_error;
typedef unsigned long devproc_selector;


#define WSDEV_WRITE_WINDOW  0
#define WSDEV_WRITE         1
#define WSDEV_HAS_DATA      2
#define WSDEV_READ          3
#define WSDEV_READ_COMPLETE 4
#define WSDEV_ACKNOWLEDGE   5

#define T_DEV_WRITE_WINDOW  WSDEV_WRITE_WINDOW
#define T_DEV_WRITE         WSDEV_WRITE
#define T_DEV_HAS_DATA      WSDEV_HAS_DATA
#define T_DEV_READ          WSDEV_READ
#define T_DEV_READ_COMPLETE WSDEV_READ_COMPLETE


#ifndef SCATTERED
#define SCATTERED 0
#undef NOT_SCATTERED
#define NOT_SCATTERED 1
#endif


typedef struct read_buf {
	unsigned short length;
	unsigned char* ptr;
} read_buf;

typedef read_buf FAR * read_bufp;
typedef read_bufp FAR * read_bufpp;

WSDMPROC( devproc_error, WSDeviceProcPtr, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));
WSDMDECL( devproc_error, WSDeviceMain, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));

typedef WSDeviceProcPtr WSDeviceUPP;
#define CallWSDeviceProc(userRoutine, thing, selector, p1, p2) (*(userRoutine))((thing), (selector), (dev_voidp)(p1), (dev_voidp)(p2))
#define NewWSDeviceProc(userRoutine) (userRoutine)

typedef WSDeviceUPP dev_main_type;
typedef dev_main_type FAR * dev_main_typep;

WS_END_EXTERN_C


#endif /* _WSDEVICE_H */




#ifndef _WSAPI_H
#define _WSAPI_H










WS_EXTERN_C

#ifndef _WSALLOC_H
#define _WSALLOC_H






#if WSINTERFACE >= 3
#if WIN64_WSTP
WSAPROC( void*, WSAllocatorProcPtr, (unsigned __int64));
#else
WSAPROC( void*, WSAllocatorProcPtr, (unsigned long));
#endif


#else

#if WIN64_WSTP
WSAPROC( WSPointer, WSAllocatorProcPtr, (unsigned __int64));
#else
WSAPROC( WSPointer, WSAllocatorProcPtr, (unsigned long));
#endif

#endif /* WSINTERFACE >= 3 */

typedef WSAllocatorProcPtr WSAllocatorUPP;
#define CallWSAllocatorProc(userRoutine, size) (*(userRoutine))((size))
#define NewWSAllocatorProc(userRoutine) (userRoutine)




#if WSINTERFACE >= 3
WSFPROC( void, WSDeallocatorProcPtr, (void*));
#else
WSFPROC( void, WSDeallocatorProcPtr, (WSPointer));
#endif /* WSINTERFACE >= 3 */

typedef WSDeallocatorProcPtr WSDeallocatorUPP;
#define CallWSDeallocatorProc(userRoutine, p) (*(userRoutine))((p))
#define NewWSDeallocatorProc(userRoutine) (userRoutine)



#endif /* _WSALLOC_H */


/* explicitly not protected by _WSALLOC_H in case WSDECL is redefined for multiple inclusion */


/* just some type-safe casts */
WSDECL( __WSProcPtr__, WSAllocatorCast,   ( WSAllocatorProcPtr f));
WSDECL( __WSProcPtr__, WSDeallocatorCast, ( WSDeallocatorProcPtr f));

WS_END_EXTERN_C



typedef WSAllocatorUPP WSAllocator;
typedef WSAllocator FAR * WSAllocatorp;
#define WSCallAllocator CallWSAllocatorProc
#define WSNewAllocator NewWSAllocatorProc

typedef WSDeallocatorUPP WSDeallocator;
typedef WSDeallocator FAR * WSDeallocatorp;
#define WSCallDeallocator CallWSDeallocatorProc
#define WSNewDeallocator NewWSDeallocatorProc

#define WSallocator WSAllocator
#define WSdeallocator WSDeallocator

#endif /* _WSAPI_H */




#ifndef _WSNTYPES_H
#define _WSNTYPES_H





#ifndef _WSNUMENV_H
#define _WSNUMENV_H








WS_EXTERN_C


#define REALBIT 4
#define REAL_MASK (1 << REALBIT)
#define XDRBIT 5
#define XDR_MASK (1 << XDRBIT)
#define BINARYBIT 7
#define BINARY_MASK (1 << BINARYBIT)
#define SIZEVARIANTBIT 6
#define SIZEVARIANT_MASK (1 << SIZEVARIANTBIT)


#define WSTK_INVALID                                          155


#define WSNE__IMPLIED_SIZE( tok, num_dispatch) ((tok) & XDR_MASK || !((tok) & SIZEVARIANT_MASK) \
		? (tok) & 0x08 ? (tok) & (0x0E + 2) : (1 << ((tok)>>1 & 0x03)) \
		: call_num_dispatch( (num_dispatch), WSNE__SIZESELECTOR((tok)), 0,0,0))

/* Range[-128, 127] */
/* 160 -> ((unsigned char)'\240') */
#define	WSTK_8BIT_SIGNED_2sCOMPLEMENT_INTEGER                 160
/* Range[0, 255] */
/* 161 -> ((unsigned char)'\241') */
#define	WSTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER               161
#define WSTK_8BIT_UNSIGNED_INTEGER WSTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER

/* Range[-32768, 32767] */
/* 162 -> ((unsigned char)'\242') */
#define	WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      162
/* Range[0, 65535] */
/* 163 -> ((unsigned char)'\243') */
#define	WSTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    163
#define	WSTK_16BIT_UNSIGNED_BIGENDIAN_INTEGER WSTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
/* 164 -> ((unsigned char)'\244') */
#define	WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      164
/* Range[0, 4294967295] */
/* 165 -> ((unsigned char)'\245') */
#define	WSTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    165
#define	WSTK_32BIT_UNSIGNED_BIGENDIAN_INTEGER WSTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
/* 166 -> ((unsigned char)'\246') */
#define	WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      166
/* Range[0, 18446744073709551615] */
/* 167 -> ((unsigned char)'\247') */
#define	WSTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    167
#define	WSTK_64BIT_UNSIGNED_BIGENDIAN_INTEGER WSTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER


/* Range[-32768, 32767] */
/* 226 -> ((unsigned char)'\342') */
#define	WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   226
/* Range[0, 65535] */
/* 227 -> ((unsigned char)'\343') */
#define	WSTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 227
#define	WSTK_16BIT_UNSIGNED_LITTLEENDIAN_INTEGER WSTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
/* 228 -> ((unsigned char)'\344') */
#define	WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   228
/* Range[0, 4294967295] */
/* 229 -> ((unsigned char)'\345') */
#define	WSTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 229
#define	WSTK_32BIT_UNSIGNED_LITTLEENDIAN_INTEGER WSTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
/* 230 -> ((unsigned char)'\346') */
#define	WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   230
/* Range[0, 18446744073709551615] */
/* 231 -> ((unsigned char)'\347') */
#define	WSTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 231
#define	WSTK_64BIT_UNSIGNED_LITTLEENDIAN_INTEGER WSTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER

/* Interval[{-3.402823e+38, 3.402823e+38}] */
/* 180 -> ((unsigned char)'\264')    10110100   */
#define	WSTK_BIGENDIAN_IEEE754_SINGLE	                      180
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
/* 182 -> ((unsigned char)'\266')    10110110   */
#define	WSTK_BIGENDIAN_IEEE754_DOUBLE	                      182

/* 184 -> ((unsigned char)'\270')    10111000   */
#define WSTK_BIGENDIAN_128BIT_DOUBLE                          184

/* Interval[{-3.402823e+38, 3.402823e+38}] */
/* 244 -> ((unsigned char)'\364')    11110100   */
#define	WSTK_LITTLEENDIAN_IEEE754_SINGLE	                  244
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
/* 246 -> ((unsigned char)'\366')    11110110   */
#define	WSTK_LITTLEENDIAN_IEEE754_DOUBLE	                  246

/* 248 -> ((unsigned char)'\370')    11111000   */
#define WSTK_LITTLEENDIAN_128BIT_DOUBLE                       248


/* Note, if the future brings...
 * #define WSTK_128BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER   ((unsigned char)'\257')
 * with  Range[0, 340282366920938463463374607431768211456 (*approximately 3.40282e+38*)]
 * the dynamic range is still a monotonically increasing function of the token value.
 * An implementation might choose to set the high varient bit to mainain this property
 * and dispatch more efficiently by avoiding overflow checks
 */

#if USE_WSNE__SELECTOR_FUNCTION
int WSNE__SELECTOR(int dtok, int stok);
#else
#define WSNE__SELECTOR( dtok, stok) \
	(((dtok) << 8) | (stok)) /* maybe should mask of high word and cast stok */
#endif
#define WSNE__SIZESELECTOR( tok) WSNE__SELECTOR( 0, tok) 
#define WSNE__INITSELECTOR (0)
#define WSNE__TOSTRINGSELECTOR( tok) WSNE__SELECTOR( WSNE__IS_REAL(tok) ? WSTKREAL : WSTKINT, tok)
#define WSNE__FROMSTRINGSELECTOR( dtok, stok) WSNE__SELECTOR( dtok, stok) 

#define WSNE__STOK( selector) ( (selector) & 0x000000FF)
#define WSNE__DTOK( selector) ( ((selector) & 0x0000FF00)>>8)

#define WSNE__IS_BINARY( tok) ((tok) & BINARY_MASK)
#define WSNE__IS_REAL( tok) ((tok) & REAL_MASK)
#define WSNE__TEXT_TOKEN( tok) (WSNE__IS_REAL( tok) ? WSTKREAL : WSTKINT)




WS_END_EXTERN_C


#endif /* _WSNUMENV_H */




#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif


/****************  Special Token types: ****************/

/* WSTK_CSHORT_P         193
   WSTK_CINT_P           194
   WSTK_CLONG_P          195
   WSTK_CFLOAT_P         209
   WSTK_CDOUBLE_P        210
   WSTK_CLONGDOUBLE_P    211 */

#define WSTK_CSHORT_P       (( BINARY_MASK | SIZEVARIANT_MASK | 1))
#define WSTK_CINT_P         (( BINARY_MASK | SIZEVARIANT_MASK | 2))
#define WSTK_CLONG_P        (( BINARY_MASK | SIZEVARIANT_MASK | 3))
#define WSTK_CFLOAT_P       (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 1))
#define WSTK_CDOUBLE_P      (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 2))
#define WSTK_CLONGDOUBLE_P  (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 3))


#define WSTK_64BIT_LITTLEENDIAN_STRUCTURE 196
#define WSTK_64BIT_BIGENDIAN_STRUCTURE    197

/* 158 -> ((unsigned char)'\236') - used in Solaris numerics definitions */
#define WSTK_128BIT_EXTENDED 158
#define WSTK_128BIT_LONGDOUBLE 158


/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
/* 218 -> ((unsigned char)'\332') */
#define WSTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED 218

/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
/* ((unsigned char)'\330') */
#define WSTK_INTEL_80BIT_EXTENDED 216

/********************  MASTIFF  ****************************/
#define WSMASTIFF_NUMERICS_ID    "mastiff"
#define WSMASTIFF_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSMASTIFF_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSMASTIFF_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSMASTIFF_CLONGDOUBLE    WSTK_128BIT_EXTENDED
#define WSMASTIFF_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSMASTIFF_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSMASTIFF_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSMASTIFF_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSMASTIFF_WSLONGDOUBLE   WSTK_128BIT_EXTENDED

/********************  JAPANESECHIN  ****************************/
#define WSJAPANESECHIN_NUMERICS_ID    "japanesechin"
#define WSJAPANESECHIN_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSJAPANESECHIN_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSJAPANESECHIN_CLONGDOUBLE    WSTK_128BIT_EXTENDED
#define WSJAPANESECHIN_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSJAPANESECHIN_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSJAPANESECHIN_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSJAPANESECHIN_WSLONGDOUBLE   WSTK_128BIT_EXTENDED

/********************  BORZOI  ****************************/
/* The borzoi numerics environment specifically does not have WSBORZOI_CLONGDOUBLE or
WSBORZOI_WSLONGDOUBLE */

#define WSBORZOI_NUMERICS_ID          "borzoi"
#define WSBORZOI_CSHORT               WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_CINT                 WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_CLONG                WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_CSIZE_T              WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_CINT64               WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBORZOI_CFLOAT               WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBORZOI_CDOUBLE              WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBORZOI_WSSHORT              WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_WSINT                WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_WSLONG               WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_WSSIZE_T             WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORZOI_WSINT64              WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBORZOI_WSFLOAT              WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBORZOI_WSDOUBLE             WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  BRIARD  ****************************/
/* The briard numerics environment purposefully does not have WSBRIARD_CLONGDOUBLE or
WSBRIARD_WSLONGDOUBLE */

#define WSBRIARD_NUMERICS_ID          "briard"
#define WSBRIARD_CSHORT               WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_CINT                 WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_CLONG                WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_CINT64               WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_CSIZE_T              WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_CFLOAT               WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBRIARD_CDOUBLE              WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBRIARD_WSSHORT              WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_WSINT                WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_WSLONG               WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_WSINT64              WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_WSSIZE_T             WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBRIARD_WSFLOAT              WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBRIARD_WSDOUBLE             WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  KEESHOND  ****************************/
#define WSKEESHOND_NUMERICS_ID    "keeshond"
#define WSKEESHOND_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSKEESHOND_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSKEESHOND_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKEESHOND_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKEESHOND_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKEESHOND_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSKEESHOND_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSKEESHOND_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKEESHOND_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  KOMONDOR  ****************************/
#define WSKOMONDOR_NUMERICS_ID    "komondor"
#define WSKOMONDOR_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSKOMONDOR_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKOMONDOR_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKOMONDOR_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSKOMONDOR_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSKOMONDOR_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSKOMONDOR_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  NORWEGIANELKHOUND  ****************************/
#define WSNORWEGIANELKHOUND_NUMERICS_ID    "norwegianelkhound"
#define WSNORWEGIANELKHOUND_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSNORWEGIANELKHOUND_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSNORWEGIANELKHOUND_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSNORWEGIANELKHOUND_CLONGDOUBLE    WSTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#define WSNORWEGIANELKHOUND_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWEGIANELKHOUND_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSNORWEGIANELKHOUND_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSNORWEGIANELKHOUND_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSNORWEGIANELKHOUND_WSLONGDOUBLE   WSTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED

/********************  NORWICHTERRIOR  ****************************/
#define WSNORWICHTERRIOR_NUMERICS_ID    "norwichterrior"
#define WSNORWICHTERRIOR_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSNORWICHTERRIOR_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSNORWICHTERRIOR_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSNORWICHTERRIOR_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSNORWICHTERRIOR_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSNORWICHTERRIOR_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSNORWICHTERRIOR_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  SAINTBERNARD  ****************************/
#define WSSAINTBERNARD_NUMERICS_ID    "saintbernarnd"
#define WSSAINTBERNARD_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSSAINTBERNARD_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSAINTBERNARD_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSAINTBERNARD_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSSAINTBERNARD_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSSAINTBERNARD_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAINTBERNARD_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSAINTBERNARD_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSAINTBERNARD_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  BERNESEMOUNTAINDOG  ****************************/
#define WSBERNESEMOUNTAINDOG_NUMERICS_ID    "bernesemountaindog"
#define WSBERNESEMOUNTAINDOG_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBERNESEMOUNTAINDOG_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBERNESEMOUNTAINDOG_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSBERNESEMOUNTAINDOG_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBERNESEMOUNTAINDOG_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBERNESEMOUNTAINDOG_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBERNESEMOUNTAINDOG_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  SETTER  ****************************/
#define WSSETTER_NUMERICS_ID    "setter"
#define WSSETTER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSSETTER_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSETTER_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSETTER_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSETTER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSSETTER_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSETTER_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSETTER_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSETTER_WSLONGDOUBLE   WSTK_INTEL_80BIT_EXTENDED

/********************  FRENCH_BULLDOG  ****************************/
#define WSFRENCH_BULLDOG_NUMERICS_ID    "french_bulldog"
#define WSFRENCH_BULLDOG_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSFRENCH_BULLDOG_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSFRENCH_BULLDOG_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSFRENCH_BULLDOG_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSFRENCH_BULLDOG_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSFRENCH_BULLDOG_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSFRENCH_BULLDOG_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSFRENCH_BULLDOG_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSFRENCH_BULLDOG_WSLONGDOUBLE   WSTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  BICHON_FRISE  ****************************/
#define WSBICHON_FRISE_NUMERICS_ID    "bichon_frise"
#define WSBICHON_FRISE_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBICHON_FRISE_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBICHON_FRISE_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBICHON_FRISE_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBICHON_FRISE_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBICHON_FRISE_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBICHON_FRISE_WSLONGDOUBLE   WSTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  HELEN  ****************************/
#define WSHELEN_NUMERICS_ID    "helen"
#define WSHELEN_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSHELEN_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSHELEN_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSHELEN_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSHELEN_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSHELEN_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSHELEN_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSHELEN_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSHELEN_WSLONGDOUBLE   WSTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  BEAGLE  ****************************/
#define WSBEAGLE_NUMERICS_ID    "beagle"
#define WSBEAGLE_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBEAGLE_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBEAGLE_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBEAGLE_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBEAGLE_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBEAGLE_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBEAGLE_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBEAGLE_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBEAGLE_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  BULLTERRIER  ****************************/
#define WSBULLTERRIER_NUMERICS_ID    "bullterrier"
#define WSBULLTERRIER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBULLTERRIER_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBULLTERRIER_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBULLTERRIER_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSBULLTERRIER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSBULLTERRIER_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBULLTERRIER_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBULLTERRIER_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBULLTERRIER_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  BORDERTERRIER  ****************************/
#define WSBORDERTERRIER_NUMERICS_ID    "borderterrier"
#define WSBORDERTERRIER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBORDERTERRIER_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBORDERTERRIER_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSBORDERTERRIER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSBORDERTERRIER_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSBORDERTERRIER_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSBORDERTERRIER_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  BASENJI  ****************************/
#define WSBASENJI_NUMERICS_ID    "basenji"
#define WSBASENJI_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSBASENJI_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBASENJI_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBASENJI_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSBASENJI_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSBASENJI_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBASENJI_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBASENJI_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBASENJI_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  SHARPEI  ****************************/
#define WSSHARPEI_NUMERICS_ID    "sharpei"
#define WSSHARPEI_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSSHARPEI_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSHARPEI_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSHARPEI_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSSHARPEI_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSSHARPEI_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHARPEI_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSHARPEI_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSHARPEI_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  TIBETANMASTIFF  ****************************/
#define WSTIBETANMASTIFF_NUMERICS_ID    "tibetanmastiff"
#define WSTIBETANMASTIFF_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSTIBETANMASTIFF_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSTIBETANMASTIFF_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSTIBETANMASTIFF_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSTIBETANMASTIFF_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSTIBETANMASTIFF_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSTIBETANMASTIFF_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSTIBETANMASTIFF_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSTIBETANMASTIFF_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  GREATDANE  ****************************/
#define WSGREATDANE_NUMERICS_ID    "greatdane"
#define WSGREATDANE_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSGREATDANE_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSGREATDANE_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSGREATDANE_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSGREATDANE_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSGREATDANE_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSGREATDANE_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  REDDOG  ****************************/
#define WSREDDOG_NUMERICS_ID    "reddog"
#define WSREDDOG_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSREDDOG_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSREDDOG_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSREDDOG_CLONGDOUBLE    WSTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#define WSREDDOG_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSREDDOG_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSREDDOG_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSREDDOG_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSREDDOG_WSLONGDOUBLE   WSTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED

/********************  AUSTRALIANCATTLEDOG  ****************************/
#define WSAUSTRALIANCATTLEDOG_NUMERICS_ID    "australiancattledog"
#define WSAUSTRALIANCATTLEDOG_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSAUSTRALIANCATTLEDOG_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSAUSTRALIANCATTLEDOG_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSAUSTRALIANCATTLEDOG_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAUSTRALIANCATTLEDOG_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSAUSTRALIANCATTLEDOG_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSAUSTRALIANCATTLEDOG_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  BOXER  ****************************/
#define WSBOXER_NUMERICS_ID    "boxer"
#define WSBOXER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBOXER_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBOXER_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSBOXER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOXER_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBOXER_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBOXER_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  AKITAINU  ****************************/
#define WSAKITAINU_NUMERICS_ID    "akitainu"
#define WSAKITAINU_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSAKITAINU_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSAKITAINU_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSAKITAINU_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSAKITAINU_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSAKITAINU_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSAKITAINU_WSLONGDOUBLE   WSTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  CHIHUAHUA  ****************************/
#define WSCHIHUAHUA_NUMERICS_ID    "chihuahua"
#define WSCHIHUAHUA_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSCHIHUAHUA_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSCHIHUAHUA_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHIHUAHUA_CLONGDOUBLE    WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHIHUAHUA_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSCHIHUAHUA_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHIHUAHUA_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSCHIHUAHUA_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHIHUAHUA_WSLONGDOUBLE   WSTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  ROTTWEILER  ****************************/
#define WSROTTWEILER_NUMERICS_ID    "rottweiler"
#define WSROTTWEILER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSROTTWEILER_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSROTTWEILER_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSROTTWEILER_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSROTTWEILER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSROTTWEILER_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSROTTWEILER_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSROTTWEILER_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSROTTWEILER_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  PHARAOHHOUND  ****************************/
#define WSPHARAOHHOUND_NUMERICS_ID    "pharaohhound"
#define WSPHARAOHHOUND_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPHARAOHHOUND_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPHARAOHHOUND_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSPHARAOHHOUND_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPHARAOHHOUND_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPHARAOHHOUND_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPHARAOHHOUND_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  TROUT  ****************************/
#define WSTROUT_NUMERICS_ID    "trout"
#define WSTROUT_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSTROUT_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSTROUT_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSTROUT_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSTROUT_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSTROUT_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSTROUT_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  PUG  ****************************/
#define WSPUG_NUMERICS_ID    "pug"
#define WSPUG_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSPUG_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPUG_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPUG_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSPUG_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSPUG_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPUG_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPUG_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPUG_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  POINTER  ****************************/
#define WSPOINTER_NUMERICS_ID    "pointer"
#define WSPOINTER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSPOINTER_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPOINTER_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPOINTER_CLONGDOUBLE    WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPOINTER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSPOINTER_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSPOINTER_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSPOINTER_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSPOINTER_WSLONGDOUBLE   WSTK_BIGENDIAN_IEEE754_DOUBLE

/********************  SAMOYED  ****************************/
#define WSSAMOYED_NUMERICS_ID    "samoyed"
#define WSSAMOYED_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSAMOYED_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSAMOYED_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSSAMOYED_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSAMOYED_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSAMOYED_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSAMOYED_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  SIBERIANHUSKY  ****************************/
#define WSSIBERIANHUSKY_NUMERICS_ID    "siberianhusky"
#define WSSIBERIANHUSKY_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSSIBERIANHUSKY_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSIBERIANHUSKY_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSIBERIANHUSKY_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSSIBERIANHUSKY_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSSIBERIANHUSKY_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSSIBERIANHUSKY_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSSIBERIANHUSKY_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSSIBERIANHUSKY_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  SHIBAINU  ****************************/
#define WSSHIBAINU_NUMERICS_ID    "shibainu"
#define WSSHIBAINU_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_CFLOAT         WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSHIBAINU_CDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSHIBAINU_CLONGDOUBLE    WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSSHIBAINU_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSSHIBAINU_WSFLOAT        WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSSHIBAINU_WSDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSSHIBAINU_WSLONGDOUBLE   WSTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  NEWFOUNDLAND  ****************************/
#define WSNEWFOUNDLAND_NUMERICS_ID    "newfoundland"
#define WSNEWFOUNDLAND_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_CLONG          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_CINT64         WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSNEWFOUNDLAND_CSIZE_T        WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSNEWFOUNDLAND_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSNEWFOUNDLAND_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSNEWFOUNDLAND_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_WSLONG         WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_WSINT64        WSTK_64BIT_BIGENDIAN_STRUCTURE
#define WSNEWFOUNDLAND_WSSIZE_T       WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSNEWFOUNDLAND_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSNEWFOUNDLAND_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSNEWFOUNDLAND_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  AFFENPINSCHER  ****************************/
#define WSAFFENPINSCHER_NUMERICS_ID    "affenpinscher"
#define WSAFFENPINSCHER_CSHORT         WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_CINT           WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_CLONG          WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_CINT64         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_CSIZE_T        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_CFLOAT         WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSAFFENPINSCHER_CDOUBLE        WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSAFFENPINSCHER_CLONGDOUBLE    WSTK_BIGENDIAN_128BIT_DOUBLE
#define WSAFFENPINSCHER_WSSHORT        WSTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_WSINT          WSTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_WSLONG         WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_WSINT64        WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_WSSIZE_T       WSTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define WSAFFENPINSCHER_WSFLOAT        WSTK_BIGENDIAN_IEEE754_SINGLE
#define WSAFFENPINSCHER_WSDOUBLE       WSTK_BIGENDIAN_IEEE754_DOUBLE
#define WSAFFENPINSCHER_WSLONGDOUBLE   WSTK_BIGENDIAN_128BIT_DOUBLE

/********************  BEAUCERON  ****************************/
#define WSBEAUCERON_NUMERICS_ID        "beauceron"
#define WSBEAUCERON_CSHORT             WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_CINT               WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_CLONG              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_CINT64             WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_CSIZE_T            WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_CFLOAT             WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBEAUCERON_CDOUBLE            WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBEAUCERON_CLONGDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBEAUCERON_WSSHORT            WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_WSINT              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_WSLONG             WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_WSINT64            WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_WSSIZE_T           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBEAUCERON_WSFLOAT            WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBEAUCERON_WSDOUBLE           WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBEAUCERON_WSLONGDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE


/********************  BERGAMASCO  ****************************/
#define WSBERGAMASCO_NUMERICS_ID       "bergamasco"
#define WSBERGAMASCO_CSHORT             WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_CINT               WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_CLONG              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_CINT64             WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_CSIZE_T            WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_CFLOAT             WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBERGAMASCO_CDOUBLE            WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBERGAMASCO_CLONGDOUBLE        WSTK_LITTLEENDIAN_128BIT_DOUBLE
#define WSBERGAMASCO_WSSHORT            WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_WSINT              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_WSLONG             WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_WSINT64            WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_WSSIZE_T           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBERGAMASCO_WSFLOAT            WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBERGAMASCO_WSDOUBLE           WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBERGAMASCO_WSLONGDOUBLE       WSTK_LITTLEENDIAN_128BIT_DOUBLE


/********************  BOERBOEL  ****************************/
#define WSBOERBOEL_NUMERICS_ID       "boerboel"
#define WSBOERBOEL_CSHORT             WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_CINT               WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_CLONG              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_CINT64             WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_CSIZE_T            WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_CFLOAT             WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBOERBOEL_CDOUBLE            WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBOERBOEL_CLONGDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBOERBOEL_WSSHORT            WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_WSINT              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_WSLONG             WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_WSINT64            WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_WSSIZE_T           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSBOERBOEL_WSFLOAT            WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSBOERBOEL_WSDOUBLE           WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSBOERBOEL_WSLONGDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE


/********************  CHINOOK  ****************************/
#define WSCHINOOK_NUMERICS_ID       "chinook"
#define WSCHINOOK_CSHORT             WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_CINT               WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_CLONG              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_CINT64             WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_CSIZE_T            WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_CFLOAT             WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSCHINOOK_CDOUBLE            WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHINOOK_CLONGDOUBLE        WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHINOOK_WSSHORT            WSTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_WSINT              WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_WSLONG             WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_WSINT64            WSTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_WSSIZE_T           WSTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define WSCHINOOK_WSFLOAT            WSTK_LITTLEENDIAN_IEEE754_SINGLE
#define WSCHINOOK_WSDOUBLE           WSTK_LITTLEENDIAN_IEEE754_DOUBLE
#define WSCHINOOK_WSLONGDOUBLE       WSTK_LITTLEENDIAN_IEEE754_DOUBLE



/********************  OLD_WIN_ENV  ****************************/
#define WSOLD_WIN_ENV_NUMERICS_ID    "Sep 13 1996, 13:46:34"
#define WSOLD_WIN_ENV_CSHORT         WSTK_CSHORT_P
#define WSOLD_WIN_ENV_CINT           WSTK_CINT_P
#define WSOLD_WIN_ENV_CLONG          WSTK_CLONG_P
#define WSOLD_WIN_ENV_CINT64         WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSOLD_WIN_ENV_CSIZE_T        WSTK_CLONG_P
#define WSOLD_WIN_ENV_CFLOAT         WSTK_CFLOAT_P
#define WSOLD_WIN_ENV_CDOUBLE        WSTK_CDOUBLE_P
#define WSOLD_WIN_ENV_CLONGDOUBLE    WSTK_CLONGDOUBLE_P
#define WSOLD_WIN_ENV_WSSHORT        WSTK_CSHORT_P
#define WSOLD_WIN_ENV_WSINT          WSTK_CINT_P
#define WSOLD_WIN_ENV_WSLONG         WSTK_CLONG_P
#define WSOLD_WIN_ENV_WSINT64        WSTK_64BIT_LITTLEENDIAN_STRUCTURE
#define WSOLD_WIN_ENV_WSSIZE_T       WSTK_CLONG_P
#define WSOLD_WIN_ENV_WSFLOAT        WSTK_CFLOAT_P
#define WSOLD_WIN_ENV_WSDOUBLE       WSTK_CDOUBLE_P
#define WSOLD_WIN_ENV_WSLONGDOUBLE   WSTK_CLONGDOUBLE_P


#define WSTK_CUCHAR  WSTK_8BIT_UNSIGNED_INTEGER
#define WSTK_WSUCHAR WSTK_8BIT_UNSIGNED_INTEGER

#if UNIX_WSTP /* hueristic that works for now */
	typedef unsigned int _uint32_nt;
	typedef signed int _sint32_nt;
#else
	typedef unsigned long _uint32_nt;
	typedef signed long _sint32_nt;
#endif


#if WINDOWS_WSTP
#if WSINTERFACE > 1
#define NEW_WIN32_NUMENV 1
#endif
#endif


/* #	define WSTP_NUMERICS_ENVIRONMENT_ID "Sep 16 1996, 23:14:20" M68KMACINTOSH_MATHLINK */
/* #	define WSTP_NUMERICS_ENVIRONMENT_ID_NUMB 33 */

/* #	define WSTP_NUMERICS_ENVIRONMENT_ID "newdog" POWERMACINTOSH_MATHLINK */
/* #	define WSTP_NUMERICS_ENVIRONMENT_ID_NUMB 24 */


#if SUN_WSTP

#if __sparc || __sparc__ || sparc
#include <sys/types.h>

#if __SUNPRO_C >= 0x301 || __SUNPRO_CC >= 0x301
#if defined(_ILP32)
#define WSTP_NUMERICS_ENVIRONMENT_ID WSMASTIFF_NUMERICS_ID

#define WSTK_CSHORT       WSMASTIFF_CSHORT
#define WSTK_CINT         WSMASTIFF_CINT
#define WSTK_CLONG        WSMASTIFF_CLONG
#define WSTK_CINT64       WSMASTIFF_CINT64
#define WSTK_CSIZE_T      WSMASTIFF_CSIZE_T
#define WSTK_CFLOAT       WSMASTIFF_CFLOAT
#define WSTK_CDOUBLE      WSMASTIFF_CDOUBLE
#define WSTK_CLONGDOUBLE  WSMASTIFF_CLONGDOUBLE

#define WSTK_WSSHORT      WSMASTIFF_WSSHORT
#define WSTK_WSINT        WSMASTIFF_WSINT
#define WSTK_WSLONG       WSMASTIFF_WSLONG
#define WSTK_WSINT64      WSMASTIFF_WSINT64
#define WSTK_WSSIZE_T     WSMASTIFF_WSSIZE_T
#define WSTK_WSFLOAT      WSMASTIFF_WSFLOAT
#define WSTK_WSDOUBLE     WSMASTIFF_WSDOUBLE
#define WSTK_WSLONGDOUBLE WSMASTIFF_WSLONGDOUBLE

#elif defined(_LP64)
#define WSTP_NUMERICS_ENVIRONMENT_ID WSJAPANESECHIN_NUMERICS_ID

#define WSTK_CSHORT       WSJAPANESECHIN_CSHORT
#define WSTK_CINT         WSJAPANESECHIN_CINT
#define WSTK_CLONG        WSJAPANESECHIN_CLONG
#define WSTK_CINT64       WSJAPANESECHIN_CINT64
#define WSTK_CSIZE_T      WSJAPANESECHIN_CSIZE_T
#define WSTK_CFLOAT       WSJAPANESECHIN_CFLOAT
#define WSTK_CDOUBLE      WSJAPANESECHIN_CDOUBLE
#define WSTK_CLONGDOUBLE  WSJAPANESECHIN_CLONGDOUBLE

#define WSTK_WSSHORT      WSJAPANESECHIN_WSSHORT
#define WSTK_WSINT        WSJAPANESECHIN_WSINT
#define WSTK_WSLONG       WSJAPANESECHIN_WSLONG
#define WSTK_WSINT64      WSJAPANESECHIN_WSINT64
#define WSTK_WSSIZE_T     WSJAPANESECHIN_WSSIZE_T
#define WSTK_WSFLOAT      WSJAPANESECHIN_WSFLOAT
#define WSTK_WSDOUBLE     WSJAPANESECHIN_WSDOUBLE
#define WSTK_WSLONGDOUBLE WSJAPANESECHIN_WSLONGDOUBLE

#endif /* _ILP32 || _LP64 */
#elif defined(__GNUC__) || defined(__GNUG__)
#if defined(_ILP32)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSKEESHOND_NUMERICS_ID

#define WSTK_CSHORT        WSKEESHOND_CSHORT
#define WSTK_CINT          WSKEESHOND_CINT
#define WSTK_CLONG         WSKEESHOND_CLONG
#define WSTK_CINT64        WSKEESHOND_CINT64
#define WSTK_CSIZE_T       WSKEESHOND_CSIZE_T
#define WSTK_CFLOAT        WSKEESHOND_CFLOAT
#define WSTK_CDOUBLE       WSKEESHOND_CDOUBLE
#define WSTK_CLONGDOULBE   WSKEESHOND_CLONGDOUBLE

#define WSTK_WSSHORT       WSKEESHOND_WSSHORT
#define WSTK_WSINT         WSKEESHOND_WSINT
#define WSTK_WSLONG        WSKEESHOND_WSLONG
#define WSTK_WSINT64       WSKEESHOND_WSINT64
#define WSTK_WSSIZE_T      WSKEESHOND_WSSIZE_T
#define WSTK_WSFLOAT       WSKEESHOND_WSFLOAT
#define WSTK_WSDOUBLE      WSKEESHOND_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSKEESHOND_WSLONGDOUBLE
#elif defined(_LP64)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSKOMONDOR_NUMERICS_ID

#define WSTK_CSHORT        WSKOMONDOR_CSHORT
#define WSTK_CINT          WSKOMONDOR_CINT
#define WSTK_CLONG         WSKOMONDOR_CLONG
#define WSTK_CINT64        WSKOMONDOR_CINT64
#define WSTK_CSIZE_T       WSKOMONDOR_CSIZE_T
#define WSTK_CFLOAT        WSKOMONDOR_CFLOAT
#define WSTK_CDOUBLE       WSKOMONDOR_CDOUBLE
#define WSTK_CLONGDOULBE   WSKOMONDOR_CLONGDOUBLE

#define WSTK_WSSHORT       WSKOMONDOR_WSSHORT
#define WSTK_WSINT         WSKOMONDOR_WSINT
#define WSTK_WSLONG        WSKOMONDOR_WSLONG
#define WSTK_WSINT64       WSKOMONDOR_WSINT64
#define WSTK_WSSIZE_T      WSKOMONDOR_WSSIZE_T
#define WSTK_WSFLOAT       WSKOMONDOR_WSFLOAT
#define WSTK_WSDOUBLE      WSKOMONDOR_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSKOMONDOR_WSLONGDOUBLE

#endif /* _ILP32 || _LP64 */
#else
#if defined(_ILP32)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBORZOI_NUMERICS_ID

#define WSTK_CSHORT        WSBORZOI_CSHORT
#define WSTK_CINT          WSBORZOI_CINT
#define WSTK_CLONG         WSBORZOI_CLONG
#define WSTK_CINT64        WSBORZOI_CINT64
#define WSTK_CSIZE_T       WSBORZOI_CSIZE_T
#define WSTK_CFLOAT        WSBORZOI_CFLOAT
#define WSTK_CDOUBLE       WSBORZOI_CDOUBLE

#define WSTK_WSSHORT       WSBORZOI_WSSHORT
#define WSTK_WSINT         WSBORZOI_WSINT
#define WSTK_WSLONG        WSBORZOI_WSLONG
#define WSTK_WSINT64       WSBORZOI_WSINT64
#define WSTK_WSSIZE_T      WSBORZOI_WSSIZE_T
#define WSTK_WSFLOAT       WSBORZOI_WSFLOAT
#define WSTK_WSDOUBLE      WSBORZOI_WSDOUBLE

#elif defined(_LP64)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBRIARD_NUMERICS_ID

#define WSTK_CSHORT        WSBRIARD_CSHORT
#define WSTK_CINT          WSBRIARD_CINT
#define WSTK_CLONG         WSBRIARD_CLONG
#define WSTK_CINT64        WSBRIARD_CINT64
#define WSTK_CSIZE_T       WSBRIARD_CSIZE_T
#define WSTK_CFLOAT        WSBRIARD_CFLOAT
#define WSTK_CDOUBLE       WSBRIARD_CDOUBLE

#define WSTK_WSSHORT       WSBRIARD_WSSHORT
#define WSTK_WSINT         WSBRIARD_WSINT
#define WSTK_WSLONG        WSBRIARD_WSLONG
#define WSTK_WSINT64       WSBRIARD_WSINT64
#define WSTK_WSSIZE_T      WSBRIARD_WSSIZE_T
#define WSTK_WSFLOAT       WSBRIARD_WSFLOAT
#define WSTK_WSDOUBLE      WSBRIARD_WSDOUBLE

#endif /* _ILP32 || _LP64 */
			/* no error directive here as the user may be
			 * using a different compiler.  Some macros
			 * simply won't be available.
			 */
#endif /* __SUNPRO_C > 0x301 || __GNUC__ */

#elif __i386 || __i386__ || i386
#if __SUNPRO_C >= 0x301 || __SUNPRO_CC >= 0x301
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSNORWEGIANELKHOUND_NUMERICS_ID

#define WSTK_CSHORT        WSNORWEGIANELKHOUND_CSHORT
#define WSTK_CINT          WSNORWEGIANELKHOUND_CINT
#define WSTK_CLONG         WSNORWEGIANELKHOUND_CLONG
#define WSTK_CINT64        WSNORWEGIANELKHOUND_CINT64
#define WSTK_CSIZE_T       WSNORWEGIANELKHOUND_CSIZE_T
#define WSTK_CFLOAT        WSNORWEGIANELKHOUND_CFLOAT
#define WSTK_CDOUBLE       WSNORWEGIANELKHOUND_CDOUBLE
#define WSTK_CLONGDOULBE   WSNORWEGIANELKHOUND_CLONGDOUBLE

#define WSTK_WSSHORT       WSNORWEGIANELKHOUND_WSSHORT
#define WSTK_WSINT         WSNORWEGIANELKHOUND_WSINT
#define WSTK_WSLONG        WSNORWEGIANELKHOUND_WSLONG
#define WSTK_WSINT64       WSNORWEGIANELKHOUND_WSINT64
#define WSTK_WSSIZE_T      WSNORWEGIANELKHOUND_WSSIZE_T
#define WSTK_WSFLOAT       WSNORWEGIANELKHOUND_WSFLOAT
#define WSTK_WSDOUBLE      WSNORWEGIANELKHOUND_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSNORWEGIANELKHOUND_WSLONGDOUBLE
#else
			/* no error directive here as the user may be
			 * using a different compiler.  Some macros
			 * simply won't be available.
			 */
#endif /* __SUNPRO_C >= 0x301 || __SUNPRO_CC >= 0x301 */

#elif __x86_64 || __x86_64__ || x86_64
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSNORWICHTERRIOR_NUMERICS_ID

#define WSTK_CSHORT        WSNORWICHTERRIOR_CSHORT
#define WSTK_CINT          WSNORWICHTERRIOR_CINT
#define WSTK_CLONG         WSNORWICHTERRIOR_CLONG
#define WSTK_CINT64        WSNORWICHTERRIOR_CINT64
#define WSTK_CSIZE_T       WSNORWICHTERRIOR_CSIZE_T
#define WSTK_CFLOAT        WSNORWICHTERRIOR_CFLOAT
#define WSTK_CDOUBLE       WSNORWICHTERRIOR_CDOUBLE
#define WSTK_CLONGDOUBLE   WSNORWICHTERRIOR_CLONGDOUBLE

#define WSTK_WSSHORT       WSNORWICHTERRIOR_WSSHORT
#define WSTK_WSINT         WSNORWICHTERRIOR_WSINT
#define WSTK_WSLONG        WSNORWICHTERRIOR_WSLONG
#define WSTK_WSINT64       WSNORWICHTERRIOR_WSINT64
#define WSTK_WSSIZE_T      WSNORWICHTERRIOR_WSSIZE_T
#define WSTK_WSFLOAT       WSNORWICHTERRIOR_WSFLOAT
#define WSTK_WSDOUBLE      WSNORWICHTERRIOR_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSNORWICHTERRIOR_WSLONGDOUBLE

#elif __SVR4 || __svr4__
#include <sys/types.h>

#if defined(_ILP32)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSSAINTBERNARD_NUMERICS_ID

#define WSTK_CSHORT        WSSAINTBERNARD_CSHORT
#define WSTK_CINT          WSSAINTBERNARD_CINT
#define WSTK_CLONG         WSSAINTBERNARD_CLONG
#define WSTK_CINT64        WSSAINTBERNARD_CINT64
#define WSTK_CSIZE_T       WSSAINTBERNARD_CSIZE_T
#define WSTK_CFLOAT        WSSAINTBERNARD_CFLOAT
#define WSTK_CDOUBLE       WSSAINTBERNARD_CDOUBLE
#define WSTK_CLONGDOULBE   WSSAINTBERNARD_CLONGDOUBLE

#define WSTK_WSSHORT       WSSAINTBERNARD_WSSHORT
#define WSTK_WSINT         WSSAINTBERNARD_WSINT
#define WSTK_WSLONG        WSSAINTBERNARD_WSLONG
#define WSTK_WSINT64       WSSAINTBERNARD_WSINT64
#define WSTK_WSSIZE_T      WSSAINTBERNARD_WSSIZE_T
#define WSTK_WSFLOAT       WSSAINTBERNARD_WSFLOAT
#define WSTK_WSDOUBLE      WSSAINTBERNARD_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSAINTBERNARD_WSLONGDOUBLE

#elif defined(_LP64)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBERNESEMOUNTAINDOG_NUMERICS_ID

#define WSTK_CSHORT        WSBERNESEMOUNTAINDOG_CSHORT
#define WSTK_CINT          WSBERNESEMOUNTAINDOG_CINT
#define WSTK_CLONG         WSBERNESEMOUNTAINDOG_CLONG
#define WSTK_CINT64        WSBERNESEMOUNTAINDOG_CINT64
#define WSTK_CSIZE_T       WSBERNESEMOUNTAINDOG_CSIZE_T
#define WSTK_CFLOAT        WSBERNESEMOUNTAINDOG_CFLOAT
#define WSTK_CDOUBLE       WSBERNESEMOUNTAINDOG_CDOUBLE
#define WSTK_CLONGDOULBE   WSBERNESEMOUNTAINDOG_CLONGDOUBLE

#define WSTK_WSSHORT       WSBERNESEMOUNTAINDOG_WSSHORT
#define WSTK_WSINT         WSBERNESEMOUNTAINDOG_WSINT
#define WSTK_WSLONG        WSBERNESEMOUNTAINDOG_WSLONG
#define WSTK_WSINT64       WSBERNESEMOUNTAINDOG_WSINT64
#define WSTK_WSSIZE_T      WSBERNESEMOUNTAINDOG_WSSIZE_T
#define WSTK_WSFLOAT       WSBERNESEMOUNTAINDOG_WSFLOAT
#define WSTK_WSDOUBLE      WSBERNESEMOUNTAINDOG_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBERNESEMOUNTAINDOG_WSLONGDOUBLE

#endif /* _ILP32 || _LP64 */

#else
/*\ syntax\ error\ */\ )
#endif

#elif (WIN32_WSTP || WIN64_WSTP) && NEW_WIN32_NUMENV
#if WIN32_WSTP
#if WSINTERFACE < 3
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSSETTER_NUMERICS_ID

#define WSTK_CSHORT        WSSETTER_CSHORT
#define WSTK_CINT          WSSETTER_CINT
#define WSTK_CLONG         WSSETTER_CLONG
#define WSTK_CINT64        WSSETTER_CINT64
#define WSTK_CSIZE_T       WSSETTER_CSIZE_T
#define WSTK_CFLOAT        WSSETTER_CFLOAT
#define WSTK_CDOUBLE       WSSETTER_CDOUBLE
#define WSTK_CLONGDOULBE   WSSETTER_CLONGDOUBLE

#define WSTK_WSSHORT       WSSETTER_WSSHORT
#define WSTK_WSINT         WSSETTER_WSINT
#define WSTK_WSLONG        WSSETTER_WSLONG
#define WSTK_WSINT64       WSSETTER_WSINT64
#define WSTK_WSSIZE_T      WSSETTER_WSSIZE_T
#define WSTK_WSFLOAT       WSSETTER_WSFLOAT
#define WSTK_WSDOUBLE      WSSETTER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSETTER_WSLONGDOUBLE
#else /* WSINTERFACE >= 3 */
#if WSINTERFACE < 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSFRENCH_BULLDOG_NUMERICS_ID

#define WSTK_CSHORT        WSFRENCH_BULLDOG_CSHORT
#define WSTK_CINT          WSFRENCH_BULLDOG_CINT
#define WSTK_CLONG         WSFRENCH_BULLDOG_CLONG
#define WSTK_CINT64        WSFRENCH_BULLDOG_CINT64
#define WSTK_CSIZE_T       WSFRENCH_BULLDOG_CSIZE_T
#define WSTK_CFLOAT        WSFRENCH_BULLDOG_CFLOAT
#define WSTK_CDOUBLE       WSFRENCH_BULLDOG_CDOUBLE
#define WSTK_CLONGDOUBLE   WSFRENCH_BULLDOG_CLONGDOUBLE

#define WSTK_WSSHORT       WSFRENCH_BULLDOG_WSSHORT
#define WSTK_WSINT         WSFRENCH_BULLDOG_WSINT
#define WSTK_WSLONG        WSFRENCH_BULLDOG_WSLONG
#define WSTK_WSINT64       WSFRENCH_BULLDOG_WSINT64
#define WSTK_WSSIZE_T      WSFRENCH_BULLDOG_WSSIZE_T
#define WSTK_WSFLOAT       WSFRENCH_BULLDOG_WSFLOAT
#define WSTK_WSDOUBLE      WSFRENCH_BULLDOG_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSFRENCH_BULLDOG_WSLONGDOUBLE
#endif /* WSINTERFACE < 4 */

#if WSINTERFACE >= 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBOERBOEL_NUMERICS_ID

#define WSTK_CSHORT        WSBOERBOEL_CSHORT
#define WSTK_CINT          WSBOERBOEL_CINT
#define WSTK_CLONG         WSBOERBOEL_CLONG
#define WSTK_CINT64        WSBOERBOEL_CINT64
#define WSTK_CSIZE_T       WSBOERBOEL_CSIZE_T
#define WSTK_CFLOAT        WSBOERBOEL_CFLOAT
#define WSTK_CDOUBLE       WSBOERBOEL_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBOERBOEL_CLONGDOUBLE

#define WSTK_WSSHORT       WSBOERBOEL_WSSHORT
#define WSTK_WSINT         WSBOERBOEL_WSINT
#define WSTK_WSLONG        WSBOERBOEL_WSLONG
#define WSTK_WSINT64       WSBOERBOEL_WSINT64
#define WSTK_WSSIZE_T      WSBOERBOEL_WSSIZE_T
#define WSTK_WSFLOAT       WSBOERBOEL_WSFLOAT
#define WSTK_WSDOUBLE      WSBOERBOEL_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBOERBOEL_WSLONGDOUBLE
#endif /* WSINTERFACE > 4 */
#endif /* WSINTERFACE < 3 */
#elif WIN64_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBICHON_FRISE_NUMERICS_ID

#define WSTK_CSHORT        WSBICHON_FRISE_CSHORT
#define WSTK_CINT          WSBICHON_FRISE_CINT
#define WSTK_CLONG         WSBICHON_FRISE_CLONG
#define WSTK_CINT64        WSBICHON_FRISE_CINT64
#define WSTK_CSIZE_T       WSBICHON_FRISE_CSIZE_T
#define WSTK_CFLOAT        WSBICHON_FRISE_CFLOAT
#define WSTK_CDOUBLE       WSBICHON_FRISE_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBICHON_FRISE_CLONGDOUBLE

#define WSTK_WSSHORT       WSBICHON_FRISE_WSSHORT
#define WSTK_WSINT         WSBICHON_FRISE_WSINT
#define WSTK_WSLONG        WSBICHON_FRISE_WSLONG
#define WSTK_WSINT64       WSBICHON_FRISE_WSINT64
#define WSTK_WSSIZE_T      WSBICHON_FRISE_WSSIZE_T
#define WSTK_WSFLOAT       WSBICHON_FRISE_WSFLOAT
#define WSTK_WSDOUBLE      WSBICHON_FRISE_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBICHON_FRISE_WSLONGDOUBLE
#endif /* WIN32_WSTP || WIN64_WSTP */

#elif ALPHA_WIN32_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSHELEN_NUMERICS_ID

#define WSTK_CSHORT        WSHELEN_CSHORT
#define WSTK_CINT          WSHELEN_CINT
#define WSTK_CLONG         WSHELEN_CLONG
#define WSTK_CINT64        WSHELEN_CINT64
#define WSTK_CSIZE_T       WSHELEN_CSIZE_T
#define WSTK_CFLOAT        WSHELEN_CFLOAT
#define WSTK_CDOUBLE       WSHELEN_CDOUBLE
#define WSTK_CLONGDOUBLE   WSHELEN_CLONGDOUBLE

#define WSTK_WSSHORT       WSHELEN_WSSHORT
#define WSTK_WSINT         WSHELEN_WSINT
#define WSTK_WSLONG        WSHELEN_WSLONG
#define WSTK_WSINT64       WSHELEN_WSINT64
#define WSTK_WSSIZE_T      WSHELEN_WSSIZE_T
#define WSTK_WSFLOAT       WSHELEN_WSFLOAT
#define WSTK_WSDOUBLE      WSHELEN_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSHELEN_WSLONGDOUBLE

#elif DARWIN_WSTP
#if PPC_DARWIN_WSTP
		/* We must use a different numerics env if we are built with gcc 3.3 or earlier. */
#if GCC_WSTP_VERSION <= 30300
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBEAGLE_NUMERICS_ID

#define WSTK_CSHORT        WSBEAGLE_CSHORT
#define WSTK_CINT          WSBEAGLE_CINT
#define WSTK_CLONG         WSBEAGLE_CLONG
#define WSTK_CINT64        WSBEAGLE_CINT64
#define WSTK_CSIZE_T       WSBEAGLE_CSIZE_T
#define WSTK_CFLOAT        WSBEAGLE_CFLOAT
#define WSTK_CDOUBLE       WSBEAGLE_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBEAGLE_CLONGDOUBLE

#define WSTK_WSSHORT       WSBEAGLE_WSSHORT
#define WSTK_WSINT         WSBEAGLE_WSINT
#define WSTK_WSLONG        WSBEAGLE_WSLONG
#define WSTK_WSINT64       WSBEAGLE_WSINT64
#define WSTK_WSSIZE_T      WSBEAGLE_WSSIZE_T
#define WSTK_WSFLOAT       WSBEAGLE_WSFLOAT
#define WSTK_WSDOUBLE      WSBEAGLE_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBEAGLE_WSLONGDOUBLE

#else
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBULLTERRIER_NUMERICS_ID

#define WSTK_CSHORT        WSBULLTERRIER_CSHORT
#define WSTK_CINT          WSBULLTERRIER_CINT
#define WSTK_CLONG         WSBULLTERRIER_CLONG
#define WSTK_CINT64        WSBULLTERRIER_CINT64
#define WSTK_CSIZE_T       WSBULLTERRIER_CSIZE_T
#define WSTK_CFLOAT        WSBULLTERRIER_CFLOAT
#define WSTK_CDOUBLE       WSBULLTERRIER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBULLTERRIER_CLONGDOUBLE

#define WSTK_WSSHORT       WSBULLTERRIER_WSSHORT
#define WSTK_WSINT         WSBULLTERRIER_WSINT
#define WSTK_WSLONG        WSBULLTERRIER_WSLONG
#define WSTK_WSINT64       WSBULLTERRIER_WSINT64
#define WSTK_WSSIZE_T      WSBULLTERRIER_WSSIZE_T
#define WSTK_WSFLOAT       WSBULLTERRIER_WSFLOAT
#define WSTK_WSDOUBLE      WSBULLTERRIER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBULLTERRIER_WSLONGDOUBLE

#endif /* GCC_WSTP_VERSION > 30300 */

#elif PPC64_DARWIN_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBORDERTERRIER_NUMERICS_ID

#define WSTK_CSHORT        WSBORDERTERRIER_CSHORT
#define WSTK_CINT          WSBORDERTERRIER_CINT
#define WSTK_CLONG         WSBORDERTERRIER_CLONG
#define WSTK_CINT64        WSBORDERTERRIER_CINT64
#define WSTK_CSIZE_T       WSBORDERTERRIER_CSIZE_T
#define WSTK_CFLOAT        WSBORDERTERRIER_CFLOAT
#define WSTK_CDOUBLE       WSBORDERTERRIER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBORDERTERRIER_CLONGDOUBLE

#define WSTK_WSSHORT       WSBORDERTERRIER_WSSHORT
#define WSTK_WSINT         WSBORDERTERRIER_WSINT
#define WSTK_WSLONG        WSBORDERTERRIER_WSLONG
#define WSTK_WSINT64       WSBORDERTERRIER_WSINT64
#define WSTK_WSSIZE_T      WSBORDERTERRIER_WSSIZE_T
#define WSTK_WSFLOAT       WSBORDERTERRIER_WSFLOAT
#define WSTK_WSDOUBLE      WSBORDERTERRIER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBORDERTERRIER_WSLONGDOUBLE

#elif X86_DARWIN_WSTP
#if WSINTERFACE < 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBASENJI_NUMERICS_ID

#define WSTK_CSHORT        WSBASENJI_CSHORT
#define WSTK_CINT          WSBASENJI_CINT
#define WSTK_CLONG         WSBASENJI_CLONG
#define WSTK_CINT64        WSBASENJI_CINT64
#define WSTK_CSIZE_T       WSBASENJI_CSIZE_T
#define WSTK_CFLOAT        WSBASENJI_CFLOAT
#define WSTK_CDOUBLE       WSBASENJI_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBASENJI_CLONGDOUBLE

#define WSTK_WSSHORT       WSBASENJI_WSSHORT
#define WSTK_WSINT         WSBASENJI_WSINT
#define WSTK_WSLONG        WSBASENJI_WSLONG
#define WSTK_WSINT64       WSBASENJI_WSINT64
#define WSTK_WSSIZE_T      WSBASENJI_WSSIZE_T
#define WSTK_WSFLOAT       WSBASENJI_WSFLOAT
#define WSTK_WSDOUBLE      WSBASENJI_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBASENJI_WSLONGDOUBLE
#endif /* WSINTERFACE < 4 */

#if WSINTERFACE >= 4
#define WSTP_NUMERICS_ENVIRONMENT_ID WSBEAUCERON_NUMERICS_ID

#define WSTK_CSHORT        WSBEAUCERON_CSHORT
#define WSTK_CINT          WSBEAUCERON_CINT
#define WSTK_CLONG         WSBEAUCERON_CLONG
#define WSTK_CINT64        WSBEAUCERON_CINT64
#define WSTK_CSIZE_T       WSBEAUCERON_CSIZE_T
#define WSTK_CFLOAT        WSBEAUCERON_CFLOAT
#define WSTK_CDOUBLE       WSBEAUCERON_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBEAUCERON_CLONGDOUBLE

#define WSTK_WSSHORT       WSBEAUCERON_WSSHORT
#define WSTK_WSINT         WSBEAUCERON_WSINT
#define WSTK_WSLONG        WSBEAUCERON_WSLONG
#define WSTK_WSINT64       WSBEAUCERON_WSINT64
#define WSTK_WSSIZE_T      WSBEAUCERON_WSSIZE_T
#define WSTK_WSFLOAT       WSBEAUCERON_WSFLOAT
#define WSTK_WSDOUBLE      WSBEAUCERON_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBEAUCERON_WSLONGDOUBLE

#endif /* WSINTERFACE >= 4 */
#elif X86_64_DARWIN_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSGREATDANE_NUMERICS_ID

#define WSTK_CSHORT        WSGREATDANE_CSHORT
#define WSTK_CINT          WSGREATDANE_CINT
#define WSTK_CLONG         WSGREATDANE_CLONG
#define WSTK_CINT64        WSGREATDANE_CINT64
#define WSTK_CSIZE_T       WSGREATDANE_CSIZE_T
#define WSTK_CFLOAT        WSGREATDANE_CFLOAT
#define WSTK_CDOUBLE       WSGREATDANE_CDOUBLE
#define WSTK_CLONGDOUBLE   WSGREATDANE_CLONGDOUBLE

#define WSTK_WSSHORT       WSGREATDANE_WSSHORT
#define WSTK_WSINT         WSGREATDANE_WSINT
#define WSTK_WSLONG        WSGREATDANE_WSLONG
#define WSTK_WSINT64       WSGREATDANE_WSINT64
#define WSTK_WSSIZE_T      WSGREATDANE_WSSIZE_T
#define WSTK_WSFLOAT       WSGREATDANE_WSFLOAT
#define WSTK_WSDOUBLE      WSGREATDANE_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSGREATDANE_WSLONGDOUBLE

#elif ARM_DARWIN_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID WSSHARPEI_NUMERICS_ID

#define WSTK_CSHORT        WSSHARPEI_CSHORT
#define WSTK_CINT          WSSHARPEI_CINT
#define WSTK_CLONG         WSSHARPEI_CLONG
#define WSTK_CINT64        WSSHARPEI_CINT64
#define WSTK_CSIZE_T       WSSHARPEI_CSIZE_T
#define WSTK_CFLOAT        WSSHARPEI_CFLOAT
#define WSTK_CDOUBLE       WSSHARPEI_CDOUBLE
#define WSTK_CLONGDOUBLE   WSSHARPEI_CLONGDOUBLE

#define WSTK_WSSHORT       WSSHARPEI_WSSHORT
#define WSTK_WSINT         WSSHARPEI_WSINT
#define WSTK_WSLONG        WSSHARPEI_WSLONG
#define WSTK_WSINT64       WSSHARPEI_WSINT64
#define WSTK_WSSIZE_T      WSSHARPEI_WSSIZE_T
#define WSTK_WSFLOAT       WSSHARPEI_WSFLOAT
#define WSTK_WSDOUBLE      WSSHARPEI_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSHARPEI_WSLONGDOUBLE
#elif ARM64_DARWIN_WSTP

#define WSTK_CSHORT        WSTIBETANMASTIFF_CSHORT
#define WSTK_CINT          WSTIBETANMASTIFF_CINT
#define WSTK_CLONG         WSTIBETANMASTIFF_CLONG
#define WSTK_CINT64        WSTIBETANMASTIFF_CINT64
#define WSTK_CSIZE_T       WSTIBETANMASTIFF_CSIZE_T
#define WSTK_CFLOAT        WSTIBETANMASTIFF_CFLOAT
#define WSTK_CDOUBLE       WSTIBETANMASTIFF_CDOUBLE
#define WSTK_CLONGDOUBLE   WSTIBETANMASTIFF_CLONGDOUBLE

#define WSTK_WSSHORT       WSTIBETANMASTIFF_WSSHORT
#define WSTK_WSINT         WSTIBETANMASTIFF_WSINT
#define WSTK_WSLONG        WSTIBETANMASTIFF_WSLONG
#define WSTK_WSINT64       WSTIBETANMASTIFF_WSINT64
#define WSTK_WSSIZE_T      WSTIBETANMASTIFF_WSSIZE_T
#define WSTK_WSFLOAT       WSTIBETANMASTIFF_WSFLOAT
#define WSTK_WSDOUBLE      WSTIBETANMASTIFF_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSTIBETANMASTIFF_WSLONGDOUBLE
#endif

#elif I86_LINUX_WSTP
#if WSINTERFACE < 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSREDDOG_NUMERICS_ID

#define WSTK_CSHORT        WSREDDOG_CSHORT
#define WSTK_CINT          WSREDDOG_CINT
#define WSTK_CLONG         WSREDDOG_CLONG
#define WSTK_CINT64        WSREDDOG_CINT64
#define WSTK_CSIZE_T       WSREDDOG_CSIZE_T
#define WSTK_CFLOAT        WSREDDOG_CFLOAT
#define WSTK_CDOUBLE       WSREDDOG_CDOUBLE
#define WSTK_CLONGDOUBLE   WSREDDOG_CLONGDOUBLE

#define WSTK_WSSHORT       WSREDDOG_WSSHORT
#define WSTK_WSINT         WSREDDOG_WSINT
#define WSTK_WSLONG        WSREDDOG_WSLONG
#define WSTK_WSINT64       WSREDDOG_WSINT64
#define WSTK_WSSIZE_T      WSREDDOG_WSSIZE_T
#define WSTK_WSFLOAT       WSREDDOG_WSFLOAT
#define WSTK_WSDOUBLE      WSREDDOG_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSREDDOG_WSLONGDOUBLE
#endif /* WSINTERFACE < 4 */

#if WSINTERFACE >= 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBERGAMASCO_NUMERICS_ID

#define WSTK_CSHORT        WSBERGAMASCO_CSHORT
#define WSTK_CINT          WSBERGAMASCO_CINT
#define WSTK_CLONG         WSBERGAMASCO_CLONG
#define WSTK_CINT64        WSBERGAMASCO_CINT64
#define WSTK_CSIZE_T       WSBERGAMASCO_CSIZE_T
#define WSTK_CFLOAT        WSBERGAMASCO_CFLOAT
#define WSTK_CDOUBLE       WSBERGAMASCO_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBERGAMASCO_CLONGDOUBLE

#define WSTK_WSSHORT       WSBERGAMASCO_WSSHORT
#define WSTK_WSINT         WSBERGAMASCO_WSINT
#define WSTK_WSLONG        WSBERGAMASCO_WSLONG
#define WSTK_WSINT64       WSBERGAMASCO_WSINT64
#define WSTK_WSSIZE_T      WSBERGAMASCO_WSSIZE_T
#define WSTK_WSFLOAT       WSBERGAMASCO_WSFLOAT
#define WSTK_WSDOUBLE      WSBERGAMASCO_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBERGAMASCO_WSLONGDOUBLE
#endif /* WSINTERFACE >= 4 */

#elif IA64_LINUX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSAUSTRALIANCATTLEDOG_NUMERICS_ID

#define WSTK_CSHORT        WSAUSTRALIANCATTLEDOG_CSHORT
#define WSTK_CINT          WSAUSTRALIANCATTLEDOG_CINT
#define WSTK_CLONG         WSAUSTRALIANCATTLEDOG_CLONG
#define WSTK_CINT64        WSAUSTRALIANCATTLEDOG_CINT64
#define WSTK_CSIZE_T       WSAUSTRALIANCATTLEDOG_CSIZE_T
#define WSTK_CFLOAT        WSAUSTRALIANCATTLEDOG_CFLOAT
#define WSTK_CDOUBLE       WSAUSTRALIANCATTLEDOG_CDOUBLE
#define WSTK_CLONGDOUBLE   WSAUSTRALIANCATTLEDOG_CLONGDOUBLE

#define WSTK_WSSHORT       WSAUSTRALIANCATTLEDOG_WSSHORT
#define WSTK_WSINT         WSAUSTRALIANCATTLEDOG_WSINT
#define WSTK_WSLONG        WSAUSTRALIANCATTLEDOG_WSLONG
#define WSTK_WSINT64       WSAUSTRALIANCATTLEDOG_WSINT64
#define WSTK_WSSIZE_T      WSAUSTRALIANCATTLEDOG_WSSIZE_T
#define WSTK_WSFLOAT       WSAUSTRALIANCATTLEDOG_WSFLOAT
#define WSTK_WSDOUBLE      WSAUSTRALIANCATTLEDOG_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSAUSTRALIANCATTLEDOG_WSLONGDOUBLE

#elif X86_64_LINUX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSBOXER_NUMERICS_ID

#define WSTK_CSHORT        WSBOXER_CSHORT
#define WSTK_CINT          WSBOXER_CINT
#define WSTK_CLONG         WSBOXER_CLONG
#define WSTK_CINT64        WSBOXER_CINT64
#define WSTK_CSIZE_T       WSBOXER_CSIZE_T
#define WSTK_CFLOAT        WSBOXER_CFLOAT
#define WSTK_CDOUBLE       WSBOXER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSBOXER_CLONGDOUBLE

#define WSTK_WSSHORT       WSBOXER_WSSHORT
#define WSTK_WSINT         WSBOXER_WSINT
#define WSTK_WSLONG        WSBOXER_WSLONG
#define WSTK_WSINT64       WSBOXER_WSINT64
#define WSTK_WSSIZE_T      WSBOXER_WSSIZE_T
#define WSTK_WSFLOAT       WSBOXER_WSFLOAT
#define WSTK_WSDOUBLE      WSBOXER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSBOXER_WSLONGDOUBLE

#elif AXP_LINUX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSAKITAINU_NUMERICS_ID

#define WSTK_CSHORT        WSAKITAINU_CSHORT
#define WSTK_CINT          WSAKITAINU_CINT
#define WSTK_CLONG         WSAKITAINU_CLONG
#define WSTK_CINT64        WSAKITAINU_CINT64
#define WSTK_CSIZE_T       WSAKITAINU_CSIZE_T
#define WSTK_CFLOAT        WSAKITAINU_CFLOAT
#define WSTK_CDOUBLE       WSAKITAINU_CDOUBLE
#define WSTK_CLONGDOUBLE   WSAKITAINU_CLONGDOUBLE

#define WSTK_WSSHORT       WSAKITAINU_WSSHORT
#define WSTK_WSINT         WSAKITAINU_WSINT
#define WSTK_WSLONG        WSAKITAINU_WSLONG
#define WSTK_WSINT64       WSAKITAINU_WSINT64
#define WSTK_WSSIZE_T      WSAKITAINU_WSSIZE_T
#define WSTK_WSFLOAT       WSAKITAINU_WSFLOAT
#define WSTK_WSDOUBLE      WSAKITAINU_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSAKITAINU_WSLONGDOUBLE

#elif ARM_LINUX_MATHLINK
#if WSINTERFACE < 4
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSCHIHUAHUA_NUMERICS_ID

#define WSTK_CSHORT        WSCHIHUAHUA_CSHORT
#define WSTK_CINT          WSCHIHUAHUA_CINT
#define WSTK_CLONG         WSCHIHUAHUA_CLONG
#define WSTK_CINT64        WSCHIHUAHUA_CINT64
#define WSTK_CSIZE_T       WSCHIHUAHUA_CSIZE_T
#define WSTK_CFLOAT        WSCHIHUAHUA_CFLOAT
#define WSTK_CDOUBLE       WSCHIHUAHUA_CDOUBLE
#define WSTK_CLONGDOUBLE   WSCHIHUAHUA_CLONGDOUBLE

#define WSTK_WSSHORT       WSCHIHUAHUA_WSSHORT
#define WSTK_WSINT         WSCHIHUAHUA_WSINT
#define WSTK_WSLONG        WSCHIHUAHUA_WSLONG
#define WSTK_WSINT64       WSCHIHUAHUA_WSINT64
#define WSTK_WSSIZE_T      WSCHIHUAHUA_WSSIZE_T
#define WSTK_WSFLOAT       WSCHIHUAHUA_WSFLOAT
#define WSTK_WSDOUBLE      WSCHIHUAHUA_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSCHIHUAHUA_WSLONGDOUBLE
#else /* WSINTERFACE >= 4 */

#define WSTP_NUMERICS_ENVIRONMENT_ID  WSCHINOOK_NUMERICS_ID

#define WSTK_CSHORT        WSCHINOOK_CSHORT
#define WSTK_CINT          WSCHINOOK_CINT
#define WSTK_CLONG         WSCHINOOK_CLONG
#define WSTK_CINT64        WSCHINOOK_CINT64
#define WSTK_CSIZE_T       WSCHINOOK_CSIZE_T
#define WSTK_CFLOAT        WSCHINOOK_CFLOAT
#define WSTK_CDOUBLE       WSCHINOOK_CDOUBLE
#define WSTK_CLONGDOUBLE   WSCHINOOK_CLONGDOUBLE

#define WSTK_WSSHORT       WSCHINOOK_WSSHORT
#define WSTK_WSINT         WSCHINOOK_WSINT
#define WSTK_WSLONG        WSCHINOOK_WSLONG
#define WSTK_WSINT64       WSCHINOOK_WSINT64
#define WSTK_WSSIZE_T      WSCHINOOK_WSSIZE_T
#define WSTK_WSFLOAT       WSCHINOOK_WSFLOAT
#define WSTK_WSDOUBLE      WSCHINOOK_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSCHINOOK_WSLONGDOUBLE
#endif /* WSINTERFACE < 4 */

#elif PPC_LINUX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSROTTWEILER_NUMERICS_ID

#define WSTK_CSHORT        WSROTTWEILER_CSHORT
#define WSTK_CINT          WSROTTWEILER_CINT
#define WSTK_CLONG         WSROTTWEILER_CLONG
#define WSTK_CINT64        WSROTTWEILER_CINT64
#define WSTK_CSIZE_T       WSROTTWEILER_CSIZE_T
#define WSTK_CFLOAT        WSROTTWEILER_CFLOAT
#define WSTK_CDOUBLE       WSROTTWEILER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSROTTWEILER_CLONGDOUBLE

#define WSTK_WSSHORT       WSROTTWEILER_WSSHORT
#define WSTK_WSINT         WSROTTWEILER_WSINT
#define WSTK_WSLONG        WSROTTWEILER_WSLONG
#define WSTK_WSINT64       WSROTTWEILER_WSINT64
#define WSTK_WSSIZE_T      WSROTTWEILER_WSSIZE_T
#define WSTK_WSFLOAT       WSROTTWEILER_WSFLOAT
#define WSTK_WSDOUBLE      WSROTTWEILER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSROTTWEILER_WSLONGDOUBLE

#elif AIX_WSTP
#if defined(__64BIT__)
#if defined(WSTP_NUMERICS_ENVIRONMENT_ID)
#undef WSTP_NUMERICS_ENVIRONMENT_ID
#endif

#ifdef __LONGDOUBLE128
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSPHARAOHHOUND_NUMERICS_ID

#define WSTK_CSHORT        WSPHARAOHHOUND_CSHORT
#define WSTK_CINT          WSPHARAOHHOUND_CINT
#define WSTK_CLONG         WSPHARAOHHOUND_CLONG
#define WSTK_CINT64        WSPHARAOHHOUND_CINT64
#define WSTK_CSIZE_T       WSPHARAOHHOUND_CSIZE_T
#define WSTK_CFLOAT        WSPHARAOHHOUND_CFLOAT
#define WSTK_CDOUBLE       WSPHARAOHHOUND_CDOUBLE
#define WSTK_CLONGDOUBLE   WSPHARAOHHOUND_CLONGDOUBLE

#define WSTK_WSSHORT       WSPHARAOHHOUND_WSSHORT
#define WSTK_WSINT         WSPHARAOHHOUND_WSINT
#define WSTK_WSLONG        WSPHARAOHHOUND_WSLONG
#define WSTK_WSINT64       WSPHARAOHHOUND_WSINT64
#define WSTK_WSSIZE_T      WSPHARAOHHOUND_WSSIZE_T
#define WSTK_WSFLOAT       WSPHARAOHHOUND_WSFLOAT
#define WSTK_WSDOUBLE      WSPHARAOHHOUND_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSPHARAOHHOUND_WSLONGDOUBLE
#else
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSTROUT_NUMERICS_ID

#define WSTK_CSHORT        WSTROUT_CSHORT
#define WSTK_CINT          WSTROUT_CINT
#define WSTK_CLONG         WSTROUT_CLONG
#define WSTK_CINT64        WSTROUT_CINT64
#define WSTK_CSIZE_T       WSTROUT_CSIZE_T
#define WSTK_CFLOAT        WSTROUT_CFLOAT
#define WSTK_CDOUBLE       WSTROUT_CDOUBLE
#define WSTK_CLONGDOUBLE   WSTROUT_CLONGDOUBLE

#define WSTK_WSSHORT       WSTROUT_WSSHORT
#define WSTK_WSINT         WSTROUT_WSINT
#define WSTK_WSLONG        WSTROUT_WSLONG
#define WSTK_WSINT64       WSTROUT_WSINT64
#define WSTK_WSSIZE_T      WSTROUT_WSSIZE_T
#define WSTK_WSFLOAT       WSTROUT_WSFLOAT
#define WSTK_WSDOUBLE      WSTROUT_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSTROUT_WSLONGDOUBLE
#endif /* __LONGDOUBLE 128 */
#else
#if defined(WSTP_NUMERICS_ENVIRONMENT_ID)
#undef WSTP_NUMERICS_ENVIRONMENT_ID
#endif

#ifdef __LONGDOUBLE128
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSPUG_NUMERICS_ID

#define WSTK_CSHORT        WSPUG_CSHORT
#define WSTK_CINT          WSPUG_CINT
#define WSTK_CLONG         WSPUG_CLONG
#define WSTK_CINT64        WSPUG_CINT64
#define WSTK_CSIZE_T       WSPUG_CSIZE_T
#define WSTK_CFLOAT        WSPUG_CFLOAT
#define WSTK_CDOUBLE       WSPUG_CDOUBLE
#define WSTK_CLONGDOUBLE   WSPUG_CLONGDOUBLE

#define WSTK_WSSHORT       WSPUG_WSSHORT
#define WSTK_WSINT         WSPUG_WSINT
#define WSTK_WSLONG        WSPUG_WSLONG
#define WSTK_WSINT64       WSPUG_WSINT64
#define WSTK_WSSIZE_T      WSPUG_WSSIZE_T
#define WSTK_WSFLOAT       WSPUG_WSFLOAT
#define WSTK_WSDOUBLE      WSPUG_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSPUG_WSLONGDOUBLE
#else
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSPOINTER_NUMERICS_ID

#define WSTK_CSHORT        WSPOINTER_CSHORT
#define WSTK_CINT          WSPOINTER_CINT
#define WSTK_CLONG         WSPOINTER_CLONG
#define WSTK_CINT64        WSPOINTER_CINT64
#define WSTK_CSIZE_T       WSPOINTER_CSIZE_T
#define WSTK_CFLOAT        WSPOINTER_CFLOAT
#define WSTK_CDOUBLE       WSPOINTER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSPOINTER_CLONGDOUBLE

#define WSTK_WSSHORT       WSPOINTER_WSSHORT
#define WSTK_WSINT         WSPOINTER_WSINT
#define WSTK_WSLONG        WSPOINTER_WSLONG
#define WSTK_WSINT64       WSPOINTER_WSINT64
#define WSTK_WSSIZE_T      WSPOINTER_WSSIZE_T
#define WSTK_WSFLOAT       WSPOINTER_WSFLOAT
#define WSTK_WSDOUBLE      WSPOINTER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSPOINTER_WSLONGDOUBLE
#endif /* __LONGDOUBLE128 */
#endif /* __64BIT__ */

#elif HPUX_WSTP
#if defined(__LP64__)
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSSAMOYED_NUMERICS_ID

#define WSTK_CSHORT        WSSAMOYED_CSHORT
#define WSTK_CINT          WSSAMOYED_CINT
#define WSTK_CLONG         WSSAMOYED_CLONG
#define WSTK_CINT64        WSSAMOYED_CINT64
#define WSTK_CSIZE_T       WSSAMOYED_CSIZE_T
#define WSTK_CFLOAT        WSSAMOYED_CFLOAT
#define WSTK_CDOUBLE       WSSAMOYED_CDOUBLE
#define WSTK_CLONGDOUBLE   WSSAMOYED_CLONGDOUBLE

#define WSTK_WSSHORT       WSSAMOYED_WSSHORT
#define WSTK_WSINT         WSSAMOYED_WSINT
#define WSTK_WSLONG        WSSAMOYED_WSLONG
#define WSTK_WSINT64       WSSAMOYED_WSINT64
#define WSTK_WSSIZE_T      WSSAMOYED_WSSIZE_T
#define WSTK_WSFLOAT       WSSAMOYED_WSFLOAT
#define WSTK_WSDOUBLE      WSSAMOYED_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSAMOYED_WSLONGDOUBLE
#else
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSSIBERIANHUSKY_NUMERICS_ID

#define WSTK_CSHORT        WSSIBERIANHUSKY_CSHORT
#define WSTK_CINT          WSSIBERIANHUSKY_CINT
#define WSTK_CLONG         WSSIBERIANHUSKY_CLONG
#define WSTK_CINT64        WSSIBERIANHUSKY_CINT64
#define WSTK_CSIZE_T       WSSIBERIANHUSKY_CSIZE_T
#define WSTK_CFLOAT        WSSIBERIANHUSKY_CFLOAT
#define WSTK_CDOUBLE       WSSIBERIANHUSKY_CDOUBLE
#define WSTK_CLONGDOUBLE   WSSIBERIANHUSKY_CLONGDOUBLE

#define WSTK_WSSHORT       WSSIBERIANHUSKY_WSSHORT
#define WSTK_WSINT         WSSIBERIANHUSKY_WSINT
#define WSTK_WSLONG        WSSIBERIANHUSKY_WSLONG
#define WSTK_WSINT64       WSSIBERIANHUSKY_WSINT64
#define WSTK_WSSIZE_T      WSSIBERIANHUSKY_WSSIZE_T
#define WSTK_WSFLOAT       WSSIBERIANHUSKY_WSFLOAT
#define WSTK_WSDOUBLE      WSSIBERIANHUSKY_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSIBERIANHUSKY_WSLONGDOUBLE
#endif /* __LP64__ */

#elif DIGITAL_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSSHIBAINU_NUMERICS_ID

#define WSTK_CSHORT        WSSHIBAINU_CSHORT
#define WSTK_CINT          WSSHIBAINU_CINT
#define WSTK_CLONG         WSSHIBAINU_CLONG
#define WSTK_CINT64        WSSHIBAINU_CINT64
#define WSTK_CSIZE_T       WSSHIBAINU_CSIZE_T
#define WSTK_CFLOAT        WSSHIBAINU_CFLOAT
#define WSTK_CDOUBLE       WSSHIBAINU_CDOUBLE
#define WSTK_CLONGDOUBLE   WSSHIBAINU_CLONGDOUBLE

#define WSTK_WSSHORT       WSSHIBAINU_WSSHORT
#define WSTK_WSINT         WSSHIBAINU_WSINT
#define WSTK_WSLONG        WSSHIBAINU_WSLONG
#define WSTK_WSINT64       WSSHIBAINU_WSINT64
#define WSTK_WSSIZE_T      WSSHIBAINU_WSSIZE_T
#define WSTK_WSFLOAT       WSSHIBAINU_WSFLOAT
#define WSTK_WSDOUBLE      WSSHIBAINU_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSSHIBAINU_WSLONGDOUBLE

#elif IRIX_WSTP
#if N32_IRIX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSNEWFOUNDLAND_NUMERICS_ID

#define WSTK_CSHORT        WSNEWFOUNDLAND_CSHORT
#define WSTK_CINT          WSNEWFOUNDLAND_CINT
#define WSTK_CLONG         WSNEWFOUNDLAND_CLONG
#define WSTK_CINT64        WSNEWFOUNDLAND_CINT64
#define WSTK_CSIZE_T       WSNEWFOUNDLAND_CSIZE_T
#define WSTK_CFLOAT        WSNEWFOUNDLAND_CFLOAT
#define WSTK_CDOUBLE       WSNEWFOUNDLAND_CDOUBLE
#define WSTK_CLONGDOUBLE   WSNEWFOUNDLAND_CLONGDOUBLE

#define WSTK_WSSHORT       WSNEWFOUNDLAND_WSSHORT
#define WSTK_WSINT         WSNEWFOUNDLAND_WSINT
#define WSTK_WSLONG        WSNEWFOUNDLAND_WSLONG
#define WSTK_WSINT64       WSNEWFOUNDLAND_WSINT64
#define WSTK_WSSIZE_T      WSNEWFOUNDLAND_WSSIZE_T
#define WSTK_WSFLOAT       WSNEWFOUNDLAND_WSFLOAT
#define WSTK_WSDOUBLE      WSNEWFOUNDLAND_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSNEWFOUNDLAND_WSLONGDOUBLE
#elif M64_IRIX_WSTP
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSAFFENPINSCHER_NUMERICS_ID

#define WSTK_CSHORT        WSAFFENPINSCHER_CSHORT
#define WSTK_CINT          WSAFFENPINSCHER_CINT
#define WSTK_CLONG         WSAFFENPINSCHER_CLONG
#define WSTK_CINT64        WSAFFENPINSCHER_CINT64
#define WSTK_CSIZE_T       WSAFFENPINSCHER_CSIZE_T
#define WSTK_CFLOAT        WSAFFENPINSCHER_CFLOAT
#define WSTK_CDOUBLE       WSAFFENPINSCHER_CDOUBLE
#define WSTK_CLONGDOUBLE   WSAFFENPINSCHER_CLONGDOUBLE

#define WSTK_WSSHORT       WSAFFENPINSCHER_WSSHORT
#define WSTK_WSINT         WSAFFENPINSCHER_WSINT
#define WSTK_WSLONG        WSAFFENPINSCHER_WSLONG
#define WSTK_WSINT64       WSAFFENPINSCHER_WSINT64
#define WSTK_WSSIZE_T      WSAFFENPINSCHER_WSSIZE_T
#define WSTK_WSFLOAT       WSAFFENPINSCHER_WSFLOAT
#define WSTK_WSDOUBLE      WSAFFENPINSCHER_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSAFFENPINSCHER_WSLONGDOUBLE
#endif

#else
#define WSTP_NUMERICS_ENVIRONMENT_ID  WSOLD_WIN_ENV_NUMERICS_ID

#define WSTK_CSHORT        WSOLD_WIN_ENV_CSHORT
#define WSTK_CINT          WSOLD_WIN_ENV_CINT
#define WSTK_CLONG         WSOLD_WIN_ENV_CLONG
#define WSTK_CINT64        WSOLD_WIN_ENV_CINT64
#define WSTK_CSIZE_T       WSOLD_WIN_ENV_CSIZE_T
#define WSTK_CFLOAT        WSOLD_WIN_ENV_CFLOAT
#define WSTK_CDOUBLE       WSOLD_WIN_ENV_CDOUBLE
#define WSTK_CLONGDOUBLE   WSOLD_WIN_ENV_CLONGDOUBLE

#define WSTK_WSSHORT       WSOLD_WIN_ENV_WSSHORT
#define WSTK_WSINT         WSOLD_WIN_ENV_WSINT
#define WSTK_WSLONG        WSOLD_WIN_ENV_WSLONG
#define WSTK_WSINT64       WSOLD_WIN_ENV_WSINT64
#define WSTK_WSSIZE_T      WSOLD_WIN_ENV_WSSIZE_T
#define WSTK_WSFLOAT       WSOLD_WIN_ENV_WSFLOAT
#define WSTK_WSDOUBLE      WSOLD_WIN_ENV_WSDOUBLE
#define WSTK_WSLONGDOUBLE  WSOLD_WIN_ENV_WSLONGDOUBLE

#endif

/* Objects of these numeric types exist in WSTP only in the numerics
 * environment and, unfortunately, in the "stack frames" of the functions that
 * put atomic numbers like WSPutInteger.  These C types are used by client
 * programs solely for type-checked access to the BinaryNumber functions.
 */
typedef unsigned char uchar_nt;
typedef uchar_nt     FAR * ucharp_nt;
typedef ucharp_nt    FAR * ucharpp_nt;

typedef short              short_nt;
typedef short_nt     FAR * shortp_nt;
typedef shortp_nt    FAR * shortpp_nt;

typedef int                int_nt;
typedef int_nt       FAR * intp_nt;
typedef intp_nt      FAR * intpp_nt;

typedef long               long_nt;
typedef long_nt      FAR * longp_nt;
typedef longp_nt     FAR * longpp_nt;

typedef wsint64            int64_nt;
typedef int64_nt         * int64p_nt;
typedef int64p_nt        * int64pp_nt;

typedef float              float_nt;
typedef float_nt     FAR * floatp_nt;
typedef floatp_nt    FAR * floatpp_nt;

typedef double             double_nt;
typedef double_nt    FAR * doublep_nt;
typedef doublep_nt   FAR * doublepp_nt;

#ifndef CC_SUPPORTS_LONG_DOUBLE
#if defined( __STDC__) || defined(__cplusplus) || ! UNIX_WSTP
#define CC_SUPPORTS_LONG_DOUBLE 1
#else

#define CC_SUPPORTS_LONG_DOUBLE WSPROTOTYPES
#endif
#endif

struct _i87extended_nt { unsigned short w[5];};

#if CC_SUPPORTS_LONG_DOUBLE
#ifndef __extended_nt__
#if WINDOWS_WSTP && (WSTK_CLONGDOUBLE != WSTK_WSLONGDOUBLE) /* subtle predicate that works for old and new windows numenvs */
#define __extended_nt__ struct _i87extended_nt
#define EXTENDED_NT_TO_I87_EXTENDED(a,b) a = b
#define I87_EXTENDED_TO_EXTENDED_NT(a,b) a = b
#else
#define __extended_nt__ long double
#define EXTENDED_NT_TO_I87_EXTENDED(a,b) \
				{ \
					int i; \
					unsigned short *c = (unsigned short *)&b; \
					for(i = 0; i < 5; i++) a.w[i] = 0; \
					for(i = 1; i < 5; i++) a.w[i] = *(c + i); \
				}
#define I87_EXTENDED_TO_EXTENDED_NT(a,b) \
				{ \
					int i; \
					unsigned short *c = (unsigned short *)&a; \
					a = 0; \
					for(i = 1; i < 5; i++) *(c + i) = b.w[i]; \
				}
#endif
#endif

	typedef __extended_nt__    wsextended_double;

	typedef __extended_nt__    extended_nt;
	typedef extended_nt  FAR * extendedp_nt;
	typedef extendedp_nt FAR * extendedpp_nt;
#endif /* CC_SUPPORTS_LONG_DOUBLE */

#endif /* _WSNTYPES_H */




#ifndef _WS0TYPES_H
#define _WS0TYPES_H




#if USING_OLD_TYPE_NAMES
typedef charp_ct ws_charp;
typedef charpp_ct ws_charpp;
typedef charppp_ct ws_charppp;
typedef ucharp_ct ws_ucharp;
typedef longp_ct ws_longp;
typedef longpp_ct ws_longpp;
typedef ulongp_ct ws_ulongp;
typedef shortp_nt ws_shortp;
typedef shortpp_nt ws_shortpp;
typedef intp_nt ws_intp;
typedef intpp_nt ws_intpp;
typedef floatp_nt ws_floatp;
typedef floatpp_nt ws_floatpp;
typedef doublep_nt ws_doublep;
typedef doublepp_nt ws_doublepp;
#if CC_SUPPORTS_LONG_DOUBLE
typedef extended_nt ws_extended;
typedef extendedp_nt ws_extendedp;
typedef extendedpp_nt ws_extendedpp;
#endif
typedef charp_ct WSBuffer;
typedef kcharp_ct WSKBuffer;
typedef charpp_ct WSBufferArray;

#endif

#endif /* _WS0TYPES_H */




WS_EXTERN_C

#ifndef _WSSTDDEV_H
#define _WSSTDDEV_H










typedef void FAR * dev_world;
typedef WSLINK dev_cookie;

typedef dev_world FAR * dev_worldp;
typedef dev_cookie FAR * dev_cookiep;


typedef  WSAllocatorUPP dev_allocator;
#define call_dev_allocator CallWSAllocatorProc
#define new_dev_allocator NewWSAllocatorProc

typedef  WSDeallocatorUPP dev_deallocator;
#define call_dev_deallocator CallWSDeallocatorProc
#define new_dev_deallocator NewWSDeallocatorProc

typedef dev_main_type world_main_type;

#define WSSTDWORLD_INIT        16
#define WSSTDWORLD_DEINIT      17
#define WSSTDWORLD_MAKE        18

#if UNIX_WSTP
#define WSSTDWORLD_GET_SIGNAL_HANDLERS      29
#define WSSTDWORLD_RELEASE_SIGNAL_HANDLERS  30
#endif

#define WSSTDWORLD_PROTOCOL        31
#define WSSTDWORLD_MODES           32
#define WSSTDWORLD_STREAMCAPACITY  33
#define WSSTDWORLD_ID              34

#define WSSTDDEV_CONNECT_READY 19
#define WSSTDDEV_CONNECT       20
#define WSSTDDEV_DESTROY       21

#define WSSTDDEV_SET_YIELDER   22
#define WSSTDDEV_GET_YIELDER   23

#define WSSTDDEV_WRITE_MSG     24
#define WSSTDDEV_HAS_MSG       25
#define WSSTDDEV_READ_MSG      26
#define WSSTDDEV_SET_HANDLER   27
#define WSSTDDEV_GET_HANDLER   28


#define T_WORLD_INIT        WSSTDWORLD_INIT
#define T_WORLD_DEINIT      WSSTDWORLD_DEINIT
#define T_WORLD_MAKE        WSSTDWORLD_MAKE
#define T_DEV_CONNECT_READY WSSTDDEV_CONNECT_READY
#define T_DEV_CONNECT       WSSTDDEV_CONNECT
#define T_DEV_DESTROY       WSSTDDEV_DESTROY

#define T_DEV_SET_YIELDER   WSSTDDEV_SET_YIELDER
#define T_DEV_GET_YIELDER   WSSTDDEV_GET_YIELDER

#define T_DEV_WRITE_MSG     WSSTDDEV_WRITE_MSG
#define T_DEV_HAS_MSG       WSSTDDEV_HAS_MSG
#define T_DEV_READ_MSG      WSSTDDEV_READ_MSG
#define T_DEV_SET_HANDLER   WSSTDDEV_SET_HANDLER
#define T_DEV_GET_HANDLER   WSSTDDEV_GET_HANDLER


typedef unsigned long dev_mode;
/* edit here and in mathlink.r */
#define NOMODE           ((dev_mode)0x0000)
#define LOOPBACKBIT      ((dev_mode)0x0001)
#define LISTENBIT        ((dev_mode)0x0002)
#define CONNECTBIT       ((dev_mode)0x0004)
#define LAUNCHBIT        ((dev_mode)0x0008)
#define PARENTCONNECTBIT ((dev_mode)0x0010)
#define READBIT          ((dev_mode)0x0020)
#define WRITEBIT         ((dev_mode)0x0040)
#define SERVERBIT        ((dev_mode)0x0080)
#define ANYMODE          (~(dev_mode)0)

typedef dev_mode FAR * dev_modep;





typedef unsigned long dev_options;

#define _DefaultOptions              ((dev_options)0x00000000)

#define _NetworkVisibleMask          ((dev_options)0x00000003)
#define _BrowseMask                  ((dev_options)0x00000010)
#define _NonBlockingMask             ((dev_options)0x00000020)
#define _InteractMask                ((dev_options)0x00000100)
#define _YieldMask                   ((dev_options)0x00000200)
#define _UseIPV6Mask                 ((dev_options)0x00010000)
#define _UseIPV4Mask                 ((dev_options)0x00020000)
#define _VersionMask                 ((dev_options)0x0F000000)
#define _UseNewTCPIPConnectionMask   ((dev_options)0x00100000)
#define _UseOldTCPIPConnectionMask   ((dev_options)0x00200000)
#define _UseUUIDTCPIPConnectionMask  ((dev_options)0x00000004)
#define _UseAnyNetworkAddressMask    ((dev_options)0x00000008)

#define _NetworkVisible              ((dev_options)0x00000000)
#define _LocallyVisible              ((dev_options)0x00000001)
#define _InternetVisible             ((dev_options)0x00000002)

#define _Browse                      ((dev_options)0x00000000)
#define _DontBrowse                  ((dev_options)0x00000010)

#define _NonBlocking                 ((dev_options)0x00000000)
#define _Blocking                    ((dev_options)0x00000020)

#define _Interact                    ((dev_options)0x00000000)
#define _DontInteract                ((dev_options)0x00000100)

#define _ForceYield                  ((dev_options)0x00000200)
#define _UseIPV6                     ((dev_options)0x00010000)
#define _UseIPV4                     ((dev_options)0x00020000)
#define _UseNewTCPIPConnection       ((dev_options)0x00100000)
#define _UseOldTCPIPConnection       ((dev_options)0x00200000)
#define _UseUUIDTCPIPConnection      ((dev_options)0x00000004)
#define _UseAnyNetworkAddress        ((dev_options)0x00000008)


#if WSINTERFACE >= 3
/* DEVICE selector and WORLD selector masks */
#define INFO_MASK (1UL << 31)
#define INFO_TYPE_MASK ((1UL << 31) - 1UL)
#define INFO_SWITCH_MASK (1UL << 30)
#define WSDEVICE_MASK INFO_MASK
#define WORLD_MASK (INFO_MASK | (1UL << 30))
#endif

/* values returned by selector WSDEVICE_TYPE */
#define UNREGISTERED_TYPE  0
#define UNIXPIPE_TYPE      1
#define UNIXSOCKET_TYPE    2
#define LOOPBACK_TYPE      5
#define WINLOCAL_TYPE      9
#define WINFMAP_TYPE       10
#define WINSHM_TYPE        11
#define SOCKET2_TYPE       12
#define GENERIC_TYPE	   13  /* Internal use only, not valid for WSDeviceInformation */
#define UNIXSHM_TYPE       14
#define INTRAPROCESS_TYPE  15

#if WSINTERFACE < 3
/* selectors */
#define WSDEVICE_TYPE 0                                       /* long */
#define WSDEVICE_NAME 1                                       /* char */
#define WSDEVICE_NAME_SIZE 2                                  /* long */
#define WSDEVICE_WORLD_ID 5                                   /* char */
#define SHM_FD                 (UNIXSHM_TYPE * 256 + 0)       /* int */
#define PIPE_FD                (UNIXPIPE_TYPE * 256 + 0)      /* int */
#define PIPE_CHILD_PID         (UNIXPIPE_TYPE * 256 + 1)      /* int */
#define SOCKET_FD              (UNIXSOCKET_TYPE * 256 + 0)    /* int */
#define INTRA_FD               (INTRAPROCESS_TYPE * 256 + 0)  /* int */
#define SOCKET_PARTNER_ADDR    (UNIXSOCKET_TYPE * 256 + 1)    /* unsigned long */
#define SOCKET_PARTNER_PORT    (UNIXSOCKET_TYPE * 256 + 2)    /* unsigned short */
#define LOOPBACK_FD            (LOOPBACK_TYPE * 256 + 2)      /* int */
#define INTRAPROCESS_FD        (INTRAPROCESS_TYPE * 256 + 0)  /* int */

#define	WINDOWS_SET_NOTIFY_WINDOW     2330 /* HWND */
#define	WINDOWS_REMOVE_NOTIFY_WINDOW  2331 /* HWND */
#define WINDOWS_READY_CONDITION       2332

/* info selectors */
#define WORLD_THISLOCATION 1        /* char */
#define WORLD_MODES 2               /* dev_mode */
#define WORLD_PROTONAME 3           /* char */
#define WORLD_STREAMCAPACITY 4      /* long */ 
#define WORLD_ID WSDEVICE_WORLD_ID    /* char */

#else /* WSINTERFACE < 3 */

/* selectors */
#define WSDEVICE_TYPE          WSDEVICE_MASK + 0UL                               /* long */
#define WSDEVICE_NAME          WSDEVICE_MASK + 1UL                               /* char */
#define WSDEVICE_NAME_SIZE     WSDEVICE_MASK + 2UL                               /* long */
#define WSDEVICE_WORLD_ID      WSDEVICE_MASK + 5UL                               /* char */
#define SHM_FD                 WSDEVICE_MASK + (UNIXSHM_TYPE * 256UL + 0UL)      /* int */
#define PIPE_FD                WSDEVICE_MASK + (UNIXPIPE_TYPE * 256UL + 0UL)     /* int */
#define PIPE_CHILD_PID         WSDEVICE_MASK + (UNIXPIPE_TYPE * 256UL + 1UL)     /* int */
#define SOCKET_FD              WSDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 0UL)   /* int */
#define INTRA_FD               WSDEVICE_MASK + (INTRAPROCESS_TYPE * 256UL + 0UL) /* int */
#define SOCKET_PARTNER_ADDR    WSDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 1UL)   /* unsigned long */
#define SOCKET_PARTNER_PORT    WSDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 2UL)   /* unsigned short */
#define LOOPBACK_FD            WSDEVICE_MASK + (LOOPBACK_TYPE * 256UL + 0UL)     /* int */
#define INTRAPROCESS_FD        WSDEVICE_MASK + (INTRAPROCESS_TYPE * 256 + 0)     /* int */

#define	WINDOWS_SET_NOTIFY_WINDOW     WSDEVICE_MASK + 2330UL /* HWND */
#define	WINDOWS_REMOVE_NOTIFY_WINDOW  WSDEVICE_MASK + 2331UL /* HWND */
#define WINDOWS_READY_CONDITION       WSDEVICE_MASK + 2332UL /* HANDLE */

/* info selectors */
#define WORLD_THISLOCATION (1UL + WORLD_MASK)        /* char */
#define WORLD_MODES (2UL + WORLD_MASK)               /* dev_mode */
#define WORLD_PROTONAME (3UL + WORLD_MASK)           /* char */
#define WORLD_STREAMCAPACITY (4UL + WORLD_MASK)      /* long */ 
#define WORLD_ID (5UL + WORLD_MASK)    /* char */
#endif /* WSINTERFACE < 3 */


#ifndef MATHLINK_DEVICE_WORLD_ID
#define MATHLINK_DEVICE_WORLD_ID (__DATE__ ", " __TIME__)
#endif


#if WSINTERFACE >= 4
#define WSDEVICE_MODE      WSDEVICE_MASK + 6UL                                    /* long */
#define WSDEVICE_OPTIONS   WSDEVICE_MASK + 7UL                                    /* long */
#endif



#define YIELDVERSION 1

typedef long devyield_result;
typedef long devyield_place;
typedef long devyield_count;
typedef unsigned long devyield_sleep;

#define INTERNAL_YIELDING 0
#define MAKE_YIELDING 1
#define CONNECT_YIELDING 2
#define READ_YIELDING 3
#define WRITE_YIELDING 4
#define DESTROY_YIELDING 5
#define READY_YIELDING 6


typedef struct WSYieldParams FAR * WSYieldParameters;


#define MAX_SLEEP (600)
typedef struct WSYieldData{
	union {long l; double d; void FAR * p;} private_data[8];
} FAR * WSYieldDataPointer;

void WSNewYieldData P(( WSYieldDataPointer ydp   /* , dev_allocator, dev_deallocator */));
void WSFreeYieldData P(( WSYieldDataPointer ydp));
WSYieldParameters WSResetYieldData P(( WSYieldDataPointer ydp, devyield_place func_id));
#if WSINTERFACE >= 3
int   WSSetYieldParameter P(( WSYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
int   WSYieldParameter P(( WSYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
#else
wsapi_result   WSSetYieldParameter P(( WSYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
wsapi_result   WSYieldParameter P(( WSYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
#endif /* WSINTERFACE >= 3 */
devyield_sleep WSSetSleepYP P(( WSYieldParameters yp, devyield_sleep sleep));
devyield_count WSSetCountYP P(( WSYieldParameters yp, devyield_count count));


enum { WSSleepParameter = 1, WSCountParameter, WSPlaceParameter};





#if WSINTERFACE >= 3
WSYPROC( int, WSYielderProcPtr, (WSLINK wslp, WSYieldParameters yp));
#else
WSYPROC( devyield_result, WSYielderProcPtr, (WSLINK wslp, WSYieldParameters yp));
#endif /* WSINTERFACE >= 3 */
typedef	WSYielderProcPtr WSDeviceYielderProcPtr;

typedef WSYielderProcPtr WSYielderUPP, WSDeviceYielderUPP;
#define NewWSYielderProc(userRoutine) (userRoutine)

#define NewWSDeviceYielderProc NewWSYielderProc

typedef  WSYielderUPP WSYieldFunctionType;

#if WSINTERFACE >= 3
typedef WSYielderUPP WSYieldFunctionObject;
#else
typedef void* WSYieldFunctionObject; /* Made change to void* for 64 bit machines */
#endif

typedef  WSYieldFunctionObject dev_yielder;
typedef dev_yielder FAR* dev_yielderp;


typedef unsigned long dev_message;
typedef dev_message FAR * dev_messagep;


#if WSINTERFACE >= 3
WSMPROC( void, WSHandlerProcPtr, (WSLINK wslp, int m, int n));
#else
WSMPROC( void, WSHandlerProcPtr, (WSLINK wslp, dev_message m, dev_message n));
#endif /* WSINTERFACE >= 3 */
typedef WSHandlerProcPtr WSDeviceHandlerProcPtr;


typedef WSHandlerProcPtr WSHandlerUPP, WSDeviceHandlerUPP;
#define NewWSHandlerProc(userRoutine) (userRoutine)

#define NewWSDeviceHandlerProc NewWSHandlerProc

typedef  WSHandlerUPP WSMessageHandlerType;

#if WSINTERFACE >= 3
typedef WSHandlerUPP WSMessageHandlerObject;
#else
#if WIN64_WSTP
typedef unsigned __int64 WSMessageHandlerObject;
#else
typedef unsigned long WSMessageHandlerObject;
#endif
#endif /* WSINTERFACE >= 3 */


typedef  WSMessageHandlerObject dev_msghandler;
typedef dev_msghandler FAR* dev_msghandlerp;



#endif /* _WSSTDDEV_H */

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif


/* explicitly not protected by _WSSTDDEV_H in case WSDECL is redefined for multiple inclusion */

/*bugcheck //should the rest of YP stuff be exported? */
WSDECL( devyield_sleep,         WSSleepYP,               ( WSYieldParameters yp));
WSDECL( devyield_count,         WSCountYP,               ( WSYieldParameters yp));

#if WSINTERFACE >= 3
WSDECL( WSYieldFunctionObject,  WSCreateYieldFunction,   ( WSEnvironment ep, WSYieldFunctionType yf, void* reserved)); /* reserved must be 0 */
#else
WSDECL( WSYieldFunctionObject,  WSCreateYieldFunction,   ( WSEnvironment ep, WSYieldFunctionType yf, WSPointer reserved)); /* reserved must be 0 */
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE > 1
#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSYieldFunctionObject,  WSCreateYieldFunction0,   ( WSEnvironment ep, WSYieldFunctionType yf, void* reserved)); /* reserved must be 0 */
#else
WSDECL( WSYieldFunctionObject,  WSCreateYieldFunction0,   ( WSEnvironment ep, WSYieldFunctionType yf, WSPointer reserved)); /* reserved must be 0 */
#endif
#endif 
#endif /* WSINTERFACE > 1 */

WSDECL( WSYieldFunctionType,    WSDestroyYieldFunction,  ( WSYieldFunctionObject yfo));

#if WSINTERFACE >= 3
WSDECL( int,        WSCallYieldFunction,     ( WSYieldFunctionObject yfo, WSLINK wslp, WSYieldParameters p));
#else
WSDECL( devyield_result,        WSCallYieldFunction,     ( WSYieldFunctionObject yfo, WSLINK wslp, WSYieldParameters p));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( WSMessageHandlerObject, WSCreateMessageHandler,  ( WSEnvironment ep, WSMessageHandlerType mh, void* reserved)); /* reserved must be 0 */
#else
WSDECL( WSMessageHandlerObject, WSCreateMessageHandler,  ( WSEnvironment ep, WSMessageHandlerType mh, WSPointer reserved)); /* reserved must be 0 */
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE > 1
#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSMessageHandlerObject, WSCreateMessageHandler0,  ( WSEnvironment ep, WSMessageHandlerType mh, void* reserved)); /* reserved must be 0 */
#else
WSDECL( WSMessageHandlerObject, WSCreateMessageHandler0,  ( WSEnvironment ep, WSMessageHandlerType mh, WSPointer reserved)); /* reserved must be 0 */
#endif
#endif
#endif /* WSINTERFACE > 1 */

WSDECL( WSMessageHandlerType,   WSDestroyMessageHandler, ( WSMessageHandlerObject mho));

#if WSINTERFACE >= 3
WSDECL( void,                   WSCallMessageHandler,    ( WSMessageHandlerObject mho, WSLINK wslp, int m, int n));
#else
WSDECL( void,                   WSCallMessageHandler,    ( WSMessageHandlerObject mho, WSLINK wslp, dev_message m, dev_message n));
#endif /* WSINTERFACE >= 3 */


/* just some type-safe casts */
WSDECL( __WSProcPtr__, WSYielderCast, ( WSYielderProcPtr yp));
WSDECL( __WSProcPtr__, WSHandlerCast, ( WSHandlerProcPtr mh));

WS_END_EXTERN_C







#ifndef WSSIGNAL_H
#define WSSIGNAL_H

#if WSINTERFACE >= 3
WSYPROC( void, WSSigHandlerProcPtr, (int signal));
#else
WSYPROC( void, WSSigHandlerProcPtr, (int_ct signal));
#endif /* WSINTERFACE >= 3 */

typedef WSSigHandlerProcPtr WSSignalHandlerType;
typedef void * WSSignalHandlerObject;

#endif /* WSSIGNAL_H */






#ifndef _WSMAKE_H
#define _WSMAKE_H












/*************** Starting WSTP ***************/

#define WSPARAMETERSIZE_R1 256
#if WSINTERFACE >= 3
#define WSPARAMETERSIZE 356
#else
#define WSPARAMETERSIZE 256
#endif

typedef char FAR * WSParametersPointer;
typedef char WSParameters[WSPARAMETERSIZE];

#define WSLoopBackOpen WSLoopbackOpen



WS_EXTERN_C
WSUPROC( void, WSUserProcPtr, (WSLINK));
WS_END_EXTERN_C

typedef WSUserProcPtr WSUserUPP;
#define NewWSUserProc(userRoutine) (userRoutine)

typedef WSUserUPP WSUserFunctionType;
typedef WSUserFunctionType FAR * WSUserFunctionTypePointer;

#if WSINTERFACE >= 3
typedef WSUserUPP WSUserFunction;
#endif /* WSINTERFACE >= 3 */


/* The following defines are
 * currently for internal use only.
 */


#define WSNetworkVisibleMask         ((unsigned long)0x00000003)          /* 00000000000000000000011 */
#define WSBrowseMask                 ((unsigned long)0x00000010)          /* 00000000000000000010000 */
#define WSNonBlockingMask            ((unsigned long)0x00000020)          /* 00000000000000000110000 */
#define WSInteractMask               ((unsigned long)0x00000100)          /* 00000000000000100000000 */
#define WSYieldMask                  ((unsigned long)0x00000200)          /* 00000000000001000000000 */
#define WSUseIPV6Mask                ((unsigned long)0x00010000)          /* 00000010000000000000000 */
#define WSUseIPV4Mask                ((unsigned long)0x00020000)          /* 00000100000000000000000 */
#define WSVersionMask                ((unsigned long)0x0000F000)          /* 00000001111000000000000 */
#define WSUseNewTCPIPConnectionMask  ((unsigned long)0x00100000)          /* 00100000000000000000000 */
#define WSUseOldTCPIPConnectionMask  ((unsigned long)0x00200000)          /* 01000000000000000000000 */
#define WSUseUUIDTCPIPConnectionMask ((unsigned long)0x00000004)          /* 00000000000000000000110 */
#define WSUseAnyNetworkAddressMask   ((unsigned long)0x00000008)          /* 00000000000000000001000 */

                                                                          /* 01100111111001100111111 */

#define WSDefaultOptions             ((unsigned long)0x00000000)
#define WSNetworkVisible             ((unsigned long)0x00000000)
#define WSLocallyVisible             ((unsigned long)0x00000001)
#define WSInternetVisible            ((unsigned long)0x00000002)

#define WSBrowse                     ((unsigned long)0x00000000)
#define WSDontBrowse                 ((unsigned long)0x00000010)

#define WSNonBlocking                ((unsigned long)0x00000000)
#define WSBlocking                   ((unsigned long)0x00000020)

#define WSInteract                   ((unsigned long)0x00000000)
#define WSDontInteract               ((unsigned long)0x00000100)

#define WSForceYield                 ((unsigned long)0x00000200)
#define WSUseIPV6                    ((unsigned long)0x00010000)
#define WSUseIPV4                    ((unsigned long)0x00020000) 

#define WSUseNewTCPIPConnection      ((unsigned long)0x00100000)
#define WSUseOldTCPIPConnection      ((unsigned long)0x00200000)
#define WSUseUUIDTCPIPConnection     ((unsigned long)0x00000004)

#define WSUseAnyNetworkAddress       ((unsigned long)0x00000008)

/* Encoding types for use with WSSetEncodingParameter */
#if WSINTERFACE >= 3
#define WSASCII_ENC		1
#define WSBYTES_ENC		2
#define WSUCS2_ENC		3
#define WSOLD_ENC		4
#define WSUTF8_ENC		5
#define WSUTF16_ENC		6
#define WSUTF32_ENC		8

#define WSTOTAL_TEXT_ENCODINGS 8
#endif

#if WSINTERFACE >= 4
#define WSLOGERROR              0
#define WSLOGWARNING            1
#define WSLOGNOTICE             2
#define WSLOGINFO               3
#define WSLOGDEBUG              4
#define WSLOGDEBUG1             5
#define WSLOGDEBUG2             6
#define WSLOGDEBUG3             7
#define WSLOGDEBUG4             8
#endif

#endif /* _WSMAKE_H */


/* explicitly not protected by _WSMAKE_H in case WSDECL is redefined for multiple inclusion */


WS_EXTERN_C

#if WSINTERFACE >= 4
WSDECL( WSEnvironmentParameter, WSNewParameters, (unsigned long rev, unsigned long apirev));
WSDECL( void,                   WSReleaseParameters, (WSEnvironmentParameter ep));
WSDECL( void, WSSetAllocParameter, (WSEnvironmentParameter ep, WSAllocator allocator, WSDeallocator deallocator));
WSDECL( long, WSSetThreadSafeLinksParameter, (WSEnvironmentParameter ep));
WSDECL( int,  WSAllocParameter,       (WSEnvironmentParameter ep, WSAllocator* allocator, WSDeallocator* deallocator));
WSDECL( long, WSSetResourceParameter, (WSEnvironmentParameter ep, const char *path));
WSDECL( long, WSSetDeviceParameter,   (WSEnvironmentParameter ep, const char *devspec));
WSDECL( long, WSErrorParameter,       (WSEnvironmentParameter ep));
WSDECL( long, WSSetEncodingParameter, (WSEnvironmentParameter ep, unsigned int etype));
WSDECL( long, WSDoNotHandleSignalParameter,    (WSEnvironmentParameter ep, int signum));
#endif /* WSINTERFACE >= 4 */

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( unsigned long, WSNewParameters,     ( char* p, unsigned long rev, unsigned long apirev));
WSDECL( void,          WSSetAllocParameter, ( char* p, WSAllocator allocator, WSDeallocator deallocator));
#else
WSDECL( ulong_ct, WSNewParameters,     ( WSParametersPointer p, ulong_ct rev, ulong_ct apirev));
WSDECL( void,     WSSetAllocParameter, ( WSParametersPointer p, WSAllocator allocator, WSDeallocator deallocator));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */


#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#if WSINTERFACE > 1 && WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( int,      WSAllocParameter,       (char* p, WSAllocator* allocator, WSDeallocator* deallocator));
WSDECL( long,     WSSetResourceParameter, (char* p, const char *path));
WSDECL( long,     WSSetDeviceParameter,   (char* p, const char *devspec));
WSDECL( long,     WSErrorParameter,       (char* p));
WSDECL( long,     WSSetEncodingParameter, (char *p, unsigned int etype));
WSDECL( long,     WSDoNotHandleSignalParameter0,    (char *p, int signum));
#else
WSDECL( long,     WSAllocParameter,       (WSParametersPointer p, WSAllocatorp allocatorp, WSDeallocatorp deallocatorp));
WSDECL( long,     WSSetResourceParameter, (WSParametersPointer p, kcharp_ct path));
WSDECL( long,     WSSetDeviceParameter,   (WSParametersPointer p, kcharp_ct devspec));
WSDECL( long,     WSErrorParameter,       (WSParametersPointer p));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE > 1 && WSINTERFACE <= 3 */

#if WSINTERFACE >= 4
WSDECL( void,          WSStopHandlingSignal, (WSEnvironment env, int signum));
WSDECL( void,          WSHandleSignal,       (WSEnvironment env, int signum));
#endif /* WSINTERFACE >= 4 */

#if WSINTERFACE == 3
WSDECL( void,          WSStopHandlingSignal0,          ( WSEnvironment env, int signum));
WSDECL( void,          WSHandleSignal0,                ( WSEnvironment env, int signum));
#endif /* WSINTERFACE == 3 */

#if WSINTERFACE >= 3
WSDECL( long,          WSSetEnvironmentData,           ( WSEnvironment env, void *cookie));
WSDECL( void *,        WSEnvironmentData,              ( WSEnvironment env));
WSDECL( int,           WSSetSignalHandler,             ( WSEnvironment env, int signum, void *so));
WSDECL( int,           WSSetSignalHandlerFromFunction, ( WSEnvironment env, int signum, WSSignalHandlerType sigfunc));
WSDECL( int,           WSUnsetSignalHandler,           ( WSEnvironment env, int signum, WSSignalHandlerType sigfunc));

WSDECL( long,          WSSetSymbolReplacement,         ( WSLINK wslp, const char *priv, int prlen, const char *pub, int pblen));
WSDECL( int,           WSClearSymbolReplacement,       ( WSLINK wslp, long index));
WSDECL( void,          WSClearAllSymbolReplacements,   ( WSLINK wslp));
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( long,          WSSetSignalHandler0, ( WSEnvironment env, int signum, WSSignalHandlerObject so));
WSDECL( long,          WSUnsetSignalHandler0, ( WSEnvironment env, int signum, WSSignalHandlerObject so));
#else
WSDECL( long,          WSSetSignalHandler0, ( WSEnvironment env, int_ct signum, WSSignalHandlerObject so));
WSDECL( long,          WSUnsetSignalHandler0, ( WSEnvironment env, int_ct signum, WSSignalHandlerObject so));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

#if WSINTERFACE >= 4
WSDECL(WSEnvironment,  WSInitialize,   ( WSEnvironmentParameter ep));
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSEnvironment, WSInitialize,   ( char* p)); /* pass in NULL */
#else
WSDECL( WSEnvironment, WSInitialize,   ( WSParametersPointer p)); /* pass in NULL */
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

WSDECL( void,          WSDeinitialize, ( WSEnvironment env));

/*************** WSTP Revsion Number/Interface Number ************/
#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( void,          WSVersionNumber0, ( WSEnvironment env, long *inumb, long *rnumb));
#else
WSDECL( void,          WSVersionNumber0, ( WSEnvironment env, longp_ct inumb, longp_ct rnumb));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

#if WSINTERFACE >= 3
WSDECL( void,          WSVersionNumbers, ( WSEnvironment env, int *inumb, int *rnumb, int *bnumb));
#endif

#if WSINTERFACE >= 4
WSDECL( int,               WSCompilerID, (WSEnvironment env, const char **id));
WSDECL( void,       WSReleaseCompilerID, (WSEnvironment env, const char *id));

WSDECL( int,           WSUCS2CompilerID, (WSEnvironment env, unsigned short **id, int *length));
WSDECL( void,   WSReleaseUCS2CompilerID, (WSEnvironment env, unsigned short *id, int length));

WSDECL( int,           WSUTF8CompilerID, (WSEnvironment env, unsigned char **id, int *length));
WSDECL( void,   WSReleaseUTF8CompilerID, (WSEnvironment env, unsigned char *id, int length));

WSDECL( int,          WSUTF16CompilerID, (WSEnvironment env, unsigned short **id, int *length));
WSDECL( void,  WSReleaseUTF16CompilerID, (WSEnvironment env, unsigned short *id, int length));

WSDECL( int,          WSUTF32CompilerID, (WSEnvironment env, unsigned int **id, int *length));
WSDECL( void,  WSReleaseUTF32CompilerID, (WSEnvironment env, unsigned int *id, int length));
#endif

/********************************************************************/

/* or, if you use WSOpenArgcArgv, ...*/

#if WSINTERFACE >= 4
WSDECL( WSEnvironment, WSBegin, (WSEnvironmentParameter ep));
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSEnvironment, WSBegin, ( char* p)); /* pass in NULL */
#else
WSDECL( WSEnvironment, WSBegin, ( WSParametersPointer p)); /* pass in NULL */
#endif /* WSINTERFACE == 3 */
#endif // WSINTERFACE <= 3

WSDECL( void,          WSEnd,   ( WSEnvironment env));

/*************** Environment Identification Interface ***************/

#if WSINTERFACE >= 3
WSDECL( int, WSSetEnvIDString, ( WSEnvironment ep, const char *environment_id)); /* APPIDSERV */
WSDECL( int, WSGetLinkedEnvIDString, (WSLINK wslp, const char **environment_id)); /* APPIDSERV */
WSDECL( void, WSReleaseEnvIDString, (WSLINK wslp, const char *environment_id));
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( long, WSSetEnvIDString0, ( WSEnvironment ep, const char *environment_id)); /* APPIDSERV */
WSDECL( long, WSGetLinkedEnvIDString0, ( WSLINK wslp, const char *environment_id)); /* APPIDSERV */
#else
WSDECL( long, WSSetEnvIDString0, ( WSEnvironment ep, kcharp_ct environment_id)); /* APPIDSERV */
WSDECL( long, WSGetLinkedEnvIDString0, ( WSLINK wslp, kcharp_ct environment_id)); /* APPIDSERV */
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

/*********************************************************************/


/**************** Network Interface List API *************************/
#if WSINTERFACE >= 3
WSDECL( char **,    WSGetNetworkAddressList, ( WSEnvironment ep, unsigned long *size ));
WSDECL( void,   WSReleaseNetworkAddressList, ( WSEnvironment ep, char **addresses, unsigned long size));

#if WSINTERFACE <= 3
/* WSDisownNetworkAddressList is for internal use only. */
WSDECL( void,   WSDisownNetworkAddressList, ( WSEnvironment ep, char **addresses, unsigned long size));
#endif

WSDECL( char **,        WSGetDomainNameList, ( WSEnvironment ep, unsigned long *size ));
WSDECL( void,       WSReleaseDomainNameList, ( WSEnvironment ep, char **dnsnames, unsigned long size));

#if WSINTERFACE <= 3
/* WSDisownDomainNameList is for internal use only */
WSDECL( void,        WSDisownDomainNameList, ( WSEnvironment ep, char **dnsnames, unsigned long size));
#endif

#endif /* WSINTERFACE >= 3 */
/*********************************************************************/


/************************* Runtime Device Inspection API ***************************/
#if WSINTERFACE >= 4
WSDECL(int, WSGetAvailableLinkProtocolNames, (WSEnvironment ep, char ***protocolNames, int *length));
WSDECL(void,     WSReleaseLinkProtocolNames, (WSEnvironment ep, char **protocolNames, int length));
#endif /* WSINTERFACE >= 4 */
/*********************************************************************/


/************************* Enumerate Open Links in an Env *************/
#if WSINTERFACE >= 4
WSDECL(int,       WSGetLinksFromEnvironment, (WSEnvironment ep, WSLINK **links, int *length));
WSDECL(void,  WSReleaseLinksFromEnvironment, (WSEnvironment ep, WSLINK *links, int length));
#endif /* WSINTERFACE >= 4 */
/*********************************************************************/



#if WSNTESTPOINTS < 1
#undef WSNTESTPOINTS
#define WSNTESTPOINTS 1
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( long, WSTestPoint1, ( WSEnvironment ep, unsigned long selector, void *p1, void *p2, long *np));
#else
WSDECL( long_et, WSTestPoint1, ( WSEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif
#if WSINTERFACE > 1

#if WSNTESTPOINTS < 2
#undef WSNTESTPOINTS
#define WSNTESTPOINTS 2
#endif

#if WSINTERFACE <= 3
WSDECL( void,    WSTestPoint2,     ( WSLINK wslp));
#endif

#if WSNTESTPOINTS < 3
#undef WSNTESTPOINTS
#define WSNTESTPOINTS 3
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( unsigned long,    WSTestPoint3,     ( WSLINK wslp));
#else
WSDECL( ulong_ct,    WSTestPoint3,     ( WSLINK wslp));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

#if WSNTESTPOINTS < 4
#undef WSNTESTPOINTS
#define WSNTESTPOINTS 4
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( unsigned long,    WSTestPoint4,     ( WSLINK wslp));
#else
WSDECL( ulong_ct,    WSTestPoint4,     ( WSLINK wslp));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE <= 3 */

#if WSINTERFACE >= 4
WSDECL( long, WSNumericsQuery, ( WSEnvironment ep, unsigned long selector, void *p1, void *p2, long *np));
#else
#if WSINTERFACE == 3
WSDECL( long, WSNumberControl0, ( WSEnvironment ep, unsigned long selector, void *p1, void *p2, long *np));
#else
WSDECL( long_et, WSNumberControl0, ( WSEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np));
#endif
#endif

#else
#if WSINTERFACE >= 3
extern long WSNumberControl0( WSEnvironment ep, unsigned long selector, void *p1, void *p2, long *np);
#else
extern long_et WSNumberControl0( WSEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np);
#endif /* WSINTERFACE >= 3 */
#endif /* WSINTERFACE > 1 */


/*************** Connection interface ***************/




#if WSINTERFACE >= 4
WSDECL( int,                            WSValid, ( WSLINK wslp));
WSDECL( WSLINK, WSCreateLinkWithExternalProtocol, ( WSEnvironment ep, dev_type dev, dev_main_type dev_main, int *errp));
#endif

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSLINK,         WSCreate0,       ( WSEnvironment ep, dev_type dev, dev_main_type dev_main, int *errp));
#else
WSDECL( WSLINK,         WSCreate0,       ( WSEnvironment ep, dev_type dev, dev_main_type dev_main, longp_ct errp));
#endif
#endif /* WSINTERFACE <= 3 */


#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSLINK,         WSMake,          ( void* ep, dev_type dev, dev_main_type dev_main, int *errp));
#else
WSDECL( WSLINK,         WSMake,          ( WSPointer ep, dev_type dev, dev_main_type dev_main, longp_ct errp));
#endif /* WSINTERFACE == 3 */

WSDECL( void,          WSDestroy,       ( WSLINK wslp, dev_typep devp, dev_main_typep dev_mainp));
WSDECL( int,           WSValid0,        ( WSLINK wslp));
#endif /* WSINTERFACE <= 3 */

#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( WSLINK,         WSLoopbackOpen0, ( WSEnvironment ep, const char *features, int *errp));
#else
WSDECL( WSLINK,         WSLoopbackOpen0, ( WSEnvironment ep, kcharp_ct features, longp_ct errp));
#endif
#endif /* WSINTERFACE <= 3 */

#if WSINTERFACE >= 4
WSDECL( char **,       WSFilterArgv,   ( WSEnvironment ep, char **argv, char **argv_end));
#else
#if WSINTERFACE == 3
WSDECL( char **,       WSFilterArgv0,   ( WSEnvironment ep, char **argv, char **argv_end));
#else
WSDECL( charpp_ct,     WSFilterArgv0,   ( WSEnvironment ep, charpp_ct argv, charpp_ct argv_end));
#endif
#endif // WSINTERFACE >= 4


#if WSINTERFACE >= 3
WSDECL( long,          WSFeatureString, ( WSLINK wslp, char *buf, long buffsize));
WSDECL( WSLINK,         WSOpenArgv,      ( WSEnvironment ep, char **argv, char **argv_end, int *errp));
WSDECL( WSLINK,         WSOpenArgcArgv,  ( WSEnvironment ep, int argc, char **argv, int *errp));
WSDECL( WSLINK,         WSOpenString,    ( WSEnvironment ep, const char *command_line, int *errp));
WSDECL( WSLINK,         WSLoopbackOpen,  ( WSEnvironment ep, int *errp));
WSDECL( int,           WSStringToArgv,  ( const char *commandline, char *buf, char **argv, int len));
WSDECL( long,          WSScanString,    ( char **argv, char ***argv_end, char **commandline, char **buf));
WSDECL( long,          WSPrintArgv,     ( char *buf, char **buf_endp, char ***argvp, char **argv_end));
#else
WSDECL( long,          WSFeatureString, ( WSLINK wslp, charp_ct buf, long buffsize));
WSDECL( WSLINK,         WSOpenArgv,      ( WSEnvironment ep, charpp_ct argv, charpp_ct argv_end, longp_ct errp));
WSDECL( WSLINK,         WSOpenString,    ( WSEnvironment ep, kcharp_ct command_line, longp_ct errp));
WSDECL( WSLINK,         WSLoopbackOpen,  ( WSEnvironment ep, longp_ct errp));
WSDECL( int_ct,        WSStringToArgv,  ( kcharp_ct commandline, charp_ct buf, charpp_ct argv, int_ct len));
WSDECL( long,          WSScanString,    ( charpp_ct argv, charppp_ct argv_end, charpp_ct commandline, charpp_ct buf));
WSDECL( long,          WSPrintArgv,     ( charp_ct buf, charpp_ct buf_endp, charppp_ct argvp, charpp_ct argv_end));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( const char *,     WSErrorMessage,  ( WSLINK wslp));
WSDECL( const char *,     WSErrorString,   ( WSEnvironment env, long err));
#else
WSDECL( kcharp_ct,     WSErrorMessage,  ( WSLINK wslp));
WSDECL( kcharp_ct,     WSErrorString,   ( WSEnvironment env, long err));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL( const unsigned short *,  WSUCS2ErrorMessage,  (WSLINK wslp, int *length));
WSDECL( const unsigned char *,   WSUTF8ErrorMessage,  (WSLINK wslp, int *length));
WSDECL( const unsigned short *,  WSUTF16ErrorMessage, (WSLINK wslp, int *length));
WSDECL( const unsigned int *,    WSUTF32ErrorMessage, (WSLINK wslp, int *length));

WSDECL(void,  WSReleaseErrorMessage,      (WSLINK wslp, const char *message));
WSDECL(void,  WSReleaseUCS2ErrorMessage,  (WSLINK wslp, const unsigned short *message, int length));
WSDECL(void,  WSReleaseUTF8ErrorMessage,  (WSLINK wslp, const unsigned char *message, int length));
WSDECL(void,  WSReleaseUTF16ErrorMessage, (WSLINK wslp, const unsigned short *message, int length));
WSDECL(void,  WSReleaseUTF32ErrorMessage, (WSLINK wslp, const unsigned int *message, int length));
#endif

#if WSINTERFACE >= 3
WSDECL( WSLINK,         WSOpen,          ( int argc, char **argv));
WSDECL( WSLINK,         WSOpenInEnv,     ( WSEnvironment env, int argc, char **argv, int *errp));

#if WSINTERFACE == 3
WSDECL( WSLINK,         WSOpenS,         ( const char *command_line));
#endif

#else
WSDECL( WSLINK,         WSOpen,          ( int_ct argc, charpp_ct argv));
WSDECL( WSLINK,         WSOpenInEnv,     ( WSEnvironment env, int_ct argc, charpp_ct argv, longp_ct errp));
WSDECL( WSLINK,         WSOpenS,         ( kcharp_ct command_line));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( WSLINK,         WSDuplicateLink,   ( WSLINK parentwslp, const char *name, int *errp ));
#else
WSDECL( WSLINK,         WSDuplicateLink,   ( WSLINK parentwslp, kcharp_ct name, longp_ct errp ));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,  WSConnect,         ( WSLINK wslp));
WSDECL( int,  WSActivate,        ( WSLINK wslp));
#else
WSDECL( wsapi_result,  WSConnect,         ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE < 3
#define WSActivate WSConnect
#endif

#ifndef __feature_setp__
#define __feature_setp__
typedef struct feature_set* feature_setp;
#endif
#if WSINTERFACE >= 3
WSDECL( int,  WSEstablish,       ( WSLINK wslp, feature_setp features));
#else
WSDECL( wsapi_result,  WSEstablish,       ( WSLINK wslp, feature_setp features));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,  WSEstablishString, ( WSLINK wslp, const char *features));
#else
WSDECL( wsapi_result,  WSEstablishString, ( WSLINK wslp, kcharp_ct features));
#endif /* WSINTERFACE >= 3 */

WSDECL( void,          WSClose,           ( WSLINK wslp));

#if WSINTERFACE >= 3
WSDECL( void,          WSSetUserData,   ( WSLINK wslp, void* data, WSUserFunction f));
WSDECL( void*,         WSUserData,      ( WSLINK wslp, WSUserFunctionType *fp));
WSDECL( void,          WSSetUserBlock,  ( WSLINK wslp, void* userblock));
WSDECL( void*,         WSUserBlock,     ( WSLINK wslp));
#else
WSDECL( void,          WSSetUserData,   ( WSLINK wslp, WSPointer data, WSUserFunctionType f));
WSDECL( WSPointer,     WSUserData,      ( WSLINK wslp, WSUserFunctionTypePointer fp));
WSDECL( void,          WSSetUserBlock,  ( WSLINK wslp, WSPointer userblock));
WSDECL( WSPointer,     WSUserBlock,     ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

/* just a type-safe cast */
WSDECL( __WSProcPtr__, WSUserCast, ( WSUserProcPtr f));


#if WSINTERFACE >= 4
WSDECL(int,            WSLogStreamToFile, (WSLINK wslp, const char *name));
WSDECL(int,       WSDisableLoggingStream, (WSLINK wslp));
WSDECL(int,        WSEnableLoggingStream, (WSLINK wslp));
WSDECL(int,    WSStopLoggingStreamToFile, (WSLINK wslp, const char *name));
WSDECL(int,          WSStopLoggingStream, (WSLINK wslp));

WSDECL(int,         WSLogFileNameForLink, (WSLINK wslp, const char **name));
WSDECL(void, WSReleaseLogFileNameForLink, (WSLINK wslp, const char *name));
#endif

/* WSLinkName returns a pointer to the link's name.
 * Links are generally named when they are created
 * and are based on information that is potentially
 * useful and is available at that time.
 * Do not attempt to deallocate the name's storage
 * through this pointer.  The storage should be
 * considered in read-only memory.
 */
#if WSINTERFACE >= 3

WSDECL( const char *, WSName,    ( WSLINK wslp));
WSDECL( const char *, WSLinkName,    ( WSLINK wslp));
#else
WSDECL( kcharp_ct, WSName,    ( WSLINK wslp));
WSDECL( kcharp_ct, WSLinkName,    ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL( const unsigned short *, WSUCS2LinkName,  (WSLINK wslp, int *length));
WSDECL( const unsigned char *,  WSUTF8LinkName,  (WSLINK wslp, int *length));
WSDECL( const unsigned short *, WSUTF16LinkName, (WSLINK wslp, int *length));
WSDECL( const unsigned int *,   WSUTF32LinkName, (WSLINK wslp, int *length));

WSDECL(void, WSReleaseLinkName,      (WSLINK wslp, const char *name));
WSDECL(void, WSReleaseUCS2LinkName,  (WSLINK wslp, const unsigned short *name, int length));
WSDECL(void, WSReleaseUTF8LinkName,  (WSLINK wslp, const unsigned char *name, int length));
WSDECL(void, WSReleaseUTF16LinkName, (WSLINK wslp, const unsigned short *name, int length));
WSDECL(void, WSReleaseUTF32LinkName, (WSLINK wslp, const unsigned int *name, int length));
#endif

WSDECL( long,      WSNumber,  ( WSLINK wslp));
#if WSINTERFACE > 1
WSDECL( long,  WSToLinkID,  ( WSLINK wslp));
WSDECL( WSLINK, WSFromLinkID, ( WSEnvironment ep, long n));
#else
extern WSLINK WSFromLinkID( WSEnvironment ep, long n);
#endif

#if WSINTERFACE >= 3
WSDECL( char *,  WSSetName, ( WSLINK wslp, const char *name));
#else
WSDECL( charp_ct,  WSSetName, ( WSLINK wslp, kcharp_ct name));
#endif /* WSINTERFACE >= 3 */

/* The following functions are
 * currently for internal use only.
 */

#if WSINTERFACE >= 3
WSDECL( void*, WSInit,   ( WSallocator alloc, WSdeallocator dealloc, void* enclosing_environment));
WSDECL( void,  WSDeinit, ( void* env));
WSDECL( void*, WSEnclosingEnvironment, ( void* ep));
WSDECL( void*, WLinkEnvironment, ( WSLINK wslp));
#else
WSDECL( WSPointer, WSInit,   ( WSallocator alloc, WSdeallocator dealloc, WSPointer enclosing_environment));
WSDECL( void,      WSDeinit, ( WSPointer env));
WSDECL( WSPointer, WSEnclosingEnvironment, ( WSPointer ep));
WSDECL( WSPointer, WLinkEnvironment, ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE >= 4
WSDECL(void, WSEnableLinkLock,  (WSLINK wslp));
WSDECL(void, WSDisableLinkLock, (WSLINK wslp));
#endif /* WSINTERFACE >= 4 */


#if WSINTERFACE >= 4
WSDECL( WSEnvironment, WSLinkEnvironment, (WSLINK wslp));
#endif /* WSINTERFACE >= 4 */


#if WSINTERFACE >= 4
WSDECL( int,  WSIsLinkLoopback, (WSLINK wslp));
#endif

/* the following two functions are for internal use only */
WSDECL( WSYieldFunctionObject, WSDefaultYieldFunction,    ( WSEnvironment env));

#if WSINTERFACE >= 3
WSDECL( int,          WSSetDefaultYieldFunction, ( WSEnvironment env, WSYieldFunctionObject yf));
#else
WSDECL( wsapi_result,          WSSetDefaultYieldFunction, ( WSEnvironment env, WSYieldFunctionObject yf));
#endif /* WSINTERFACE >= 3 */


WS_END_EXTERN_C






#ifndef WSLINKSERVER_H
#define WSLINKSERVER_H






#if WSINTERFACE >= 4

WS_EXTERN_C

typedef void * WSLinkServer;

typedef void (*WSNewLinkCallbackFunction)(WSLinkServer server, WSLINK link);

WSDECL(WSLinkServer, WSNewLinkServer, (WSEnvironment env, void *context, int *error));

WSDECL(WSLinkServer, WSNewLinkServerWithPort, (WSEnvironment env, unsigned short port, void *context,
    int *error));

WSDECL(WSLinkServer, WSNewLinkServerWithPortAndInterface, (WSEnvironment env, unsigned short port, const char *iface,
    void *context, int *error));

WSDECL(void, WSShutdownLinkServer, (WSLinkServer server));

WSDECL(void, WSRegisterCallbackFunctionWithLinkServer, (WSLinkServer server, WSNewLinkCallbackFunction function));

WSDECL(WSLINK, WSWaitForNewLinkFromLinkServer, (WSLinkServer server, int *error));

WSDECL(unsigned short, WSPortFromLinkServer, (WSLinkServer server, int *error));

WSDECL(const char *, WSInterfaceFromLinkServer, (WSLinkServer server, int *error));

WSDECL(void *, WSContextFromLinkServer, (WSLinkServer server, int *error));

WSDECL(void, WSReleaseInterfaceFromLinkServer, (WSLinkServer server, const char *iface));

WS_END_EXTERN_C

#endif /* WSINTERFACE >= 4 */

#endif /* WSLINKSERVER_H */





#ifndef WSSERVICEDISCOVERYAPI_H
#define WSSERVICEDISCOVERYAPI_H







#if WSINTERFACE >= 4

WS_EXTERN_C

#define WSSDADDSERVICE      0x0001
#define WSSDREMOVESERVICE   0x0002
#define WSSDBROWSEERROR     0x0003
#define WSSDRESOLVEERROR    0x0004
#define WSSDREGISTERERROR   0x0005
#define WSSDMORECOMING      0x0010
#define WSSDNAMECONFLICT    0x0007

typedef void * WSServiceRef;

typedef void (*WSBrowseCallbackFunction)(WSEnvironment env, WSServiceRef ref, int flag, 
	const char *serviceName, void *context);

WSDECL(int,     WSBrowseForLinkServices, (WSEnvironment env, 
    WSBrowseCallbackFunction callbackFunction, const char *serviceProtocol, 
    const char *domain, void *context, WSServiceRef *ref));

WSDECL(void, WSStopBrowsingForLinkServices, (WSEnvironment env, WSServiceRef ref));

typedef void (*WSResolveCallbackFunction)(WSEnvironment env, WSServiceRef ref, const char *serviceName, 
	const char *linkName, const char *protocol, int options, void *context);

WSDECL(int, WSResolveLinkService, (WSEnvironment env,
    WSResolveCallbackFunction, const char *serviceProtocol, 
    const char *serviceName, void *context, WSServiceRef *ref));

WSDECL(void, WSStopResolvingLinkService, (WSEnvironment env, WSServiceRef ref));

typedef void (*WSRegisterCallbackFunction)(WSEnvironment env, WSServiceRef ref, int flag, const char *serviceName, 
	void *context);

WSDECL(WSLINK, WSRegisterLinkServiceWithPortAndHostname, (WSEnvironment env, const char *serviceProtocol,
    const char *serviceName, unsigned short port, const char *hostname, WSRegisterCallbackFunction function,
    const char *domain, void *context, WSServiceRef *ref, int *error));

WSDECL(WSLINK, WSRegisterLinkServiceWithHostname, (WSEnvironment env, const char *serviceProtocol,
    const char *serviceName, const char *hostname, WSRegisterCallbackFunction function,
    const char *domain, void *context, WSServiceRef *ref, int *error));

WSDECL(WSLINK, WSRegisterLinkService, (WSEnvironment env, const char *serviceProtocol,
    const char *serviceName, WSRegisterCallbackFunction function,
    const char *domain, void *context, WSServiceRef *, int *error));

WSDECL(WSLINK, WSRegisterLinkServiceUsingLinkProtocol, (WSEnvironment env, const char *serviceProtocol,
	const char *serviceName, unsigned short port, const char *hostname, const char *protocol, 
	WSRegisterCallbackFunction function, const char *domain, void *context, WSServiceRef *ref, int *error));

WSDECL(void, WSRegisterLinkServiceFromLinkServer, (WSEnvironment env, const char *serviceProtocol,
    const char *serviceName, WSLinkServer server, WSRegisterCallbackFunction function, const char *domain,
    void *context, WSServiceRef *ref, int *error));

WSDECL(void, WSStopRegisteringLinkService, (WSEnvironment env, WSServiceRef ref));

WSDECL(void, WSStopRegisteringLinkServiceForLink, (WSEnvironment env, WSLINK link, WSServiceRef ref));

WSDECL(const char *, WSServiceProtocolFromReference, (WSEnvironment env, WSServiceRef ref));

WS_END_EXTERN_C

#endif /* WSINTERFACE >= 4 */

#endif /* end of include guard: WSSERVICEDISCOVERYAPI_H */




#ifndef _WSERRORS_H
#define _WSERRORS_H




/*************** WSTP errors ***************/
/*
 * When some problem is detected within WSTP, routines
 * will return a simple indication of failure and store
 * an error code internally. (For routines that have nothing
 * else useful to return, success is indicated by returning
 * non-zero and failure by returning 0.)  WSError() returns
 * the current error code;  WSErrorMessage returns an English
 * language description of the error.
 * The error WSEDEAD is irrecoverable.  For the others, WSClearError()
 * will reset the error code to WSEOK.
 */



#ifndef _WSERRNO_H
#define _WSERRNO_H


#define WSEUNKNOWN          -1
#define WSEOK                0
#define WSEDEAD              1
#define WSEGBAD              2
#define WSEGSEQ              3
#define WSEPBTK              4
#define WSEPSEQ              5
#define WSEPBIG              6
#define WSEOVFL              7
#define WSEMEM               8
#define WSEACCEPT            9
#define WSECONNECT          10
#define WSECLOSED           11
#define WSEDEPTH            12  /* internal error */
#define WSENODUPFCN         13  /* stream cannot be duplicated */

#define WSENOACK            15  /* */
#define WSENODATA           16  /* */
#define WSENOTDELIVERED     17  /* */
#define WSENOMSG            18  /* */
#define WSEFAILED           19  /* */

#define WSEGETENDEXPR       20
#define WSEPUTENDPACKET     21 /* unexpected call of WSEndPacket */
                               /* currently atoms aren't
                                * counted on the way out so this error is raised only when
                                * WSEndPacket is called in the midst of an atom
                                */
#define WSENEXTPACKET       22
#define WSEUNKNOWNPACKET    23
#define WSEGETENDPACKET     24
#define WSEABORT            25
#define WSEMORE             26 /* internal error */
#define WSENEWLIB           27
#define WSEOLDLIB           28
#define WSEBADPARAM         29
#define WSENOTIMPLEMENTED   30


#define WSEINIT             32
#define WSEARGV             33
#define WSEPROTOCOL         34
#define WSEMODE             35
#define WSELAUNCH           36
#define WSELAUNCHAGAIN      37
#define WSELAUNCHSPACE      38
#define WSENOPARENT         39
#define WSENAMETAKEN        40
#define WSENOLISTEN         41
#define WSEBADNAME          42
#define WSEBADHOST          43
#define WSERESOURCE         44  /* a required resource was missing */
#define WSELAUNCHFAILED     45
#define WSELAUNCHNAME       46
#if WSINTERFACE < 3
#define WSELAST WSELAUNCHNAME /* for internal use only */
#elif WSINTERFACE == 3
#define WSEPDATABAD         47
#define WSEPSCONVERT        48
#define WSEGSCONVERT        49
#define WSENOTEXE           50
#define WSESYNCOBJECTMAKE   51
#define WSEBACKOUT          52
#define WSELAST WSEBACKOUT
#else /* WSINTERFACE >= 4 */
#define WSEPDATABAD         47
#define WSEPSCONVERT        48
#define WSEGSCONVERT        49
#define WSENOTEXE           50
#define WSESYNCOBJECTMAKE   51
#define WSEBACKOUT          52
#define WSEBADOPTSYM        53
#define WSEBADOPTSTR        54
#define WSENEEDBIGGERBUFFER 55
#define WSEBADNUMERICSID    56
#define WSESERVICENOTAVAILABLE 57
#define WSEBADARGUMENT      58
#define WSEBADDISCOVERYHOSTNAME         59
#define WSEBADDISCOVERYDOMAINNAME       60
#define WSEBADSERVICENAME               61
#define WSEBADDISCOVERYSTATE            62
#define WSEBADDISCOVERYFLAGS            63
#define WSEDISCOVERYNAMECOLLISION       64
#define WSEBADSERVICEDISCOVERY          65
#define WSELAST WSESERVICENOTAVAILABLE
#endif

#define WSETRACEON         996  /* */
#define WSETRACEOFF        997  /* */
#define WSEDEBUG           998  /* */
#define WSEASSERT          999  /* an internal assertion failed */
#define WSEUSER           1000  /* start of user defined errors */


#endif /* _WSERRNO_H */




#endif /* _WSERRORS_H */

/* explicitly not protected by _WSERRORS_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C
#if WSINTERFACE >= 3
WSDECL( int,  WSError,        ( WSLINK wslp));
WSDECL( int,  WSClearError,   ( WSLINK wslp));
WSDECL( int,  WLSetError,     ( WSLINK wslp, int err));
#else
WSDECL( wsapi_error,   WSError,        ( WSLINK wslp));
WSDECL( wsapi_result,  WSClearError,   ( WSLINK wslp));
WSDECL( wsapi_result,  WLSetError,     ( WSLINK wslp, wsapi_error err));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C




#ifndef _WSYLDMSG_H
#define _WSYLDMSG_H









enum {	WSTerminateMessage = 1, WSInterruptMessage, WSAbortMessage,
	WSEndPacketMessage, WSSynchronizeMessage, WSImDyingMessage,
	WSWaitingAcknowledgment, WSMarkTopLevelMessage, WSLinkClosingMessage,
	WSAuthenticateFailure, WSSuspendActivitiesMessage, WSResumeActivitiesMessage,
	WSFirstUserMessage = 128, WSLastUserMessage = 255 };

typedef unsigned long devinfo_selector;


#endif /* _WSYLDMSG_H */

/* explicitly not protected by _WSYLDMSG_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#if WSINTERFACE >= 3
WSDECL( int,   WSPutMessage,   ( WSLINK wslp, int  msg));
WSDECL( int,   WSGetMessage,   ( WSLINK wslp, int *mp, int *np));
WSDECL( int,   WSMessageReady, ( WSLINK wslp));
#else
WSDECL( wsapi_result,   WSGetMessage,   ( WSLINK wslp, dev_messagep mp, dev_messagep np));
WSDECL( wsapi_result,   WSPutMessage,   ( WSLINK wslp, dev_message  msg));
WSDECL( wsapi_result,   WSMessageReady, ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,   WSPutMessageWithArg, ( WSLINK wslp, int msg, int arg));
#endif


#if WSINTERFACE >= 3
WSDECL( WSMessageHandlerObject, WSGetMessageHandler,    ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */
WSDECL( WSMessageHandlerObject, WSMessageHandler,    ( WSLINK wslp));

#if WSINTERFACE >= 3
WSDECL( WSYieldFunctionObject,  WSGetYieldFunction,     ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */
WSDECL( WSYieldFunctionObject,  WSYieldFunction,     ( WSLINK wslp));

#if WSINTERFACE >= 3
WSDECL( int,  WSSetMessageHandler, ( WSLINK wslp, WSMessageHandlerObject h));
WSDECL( int,  WSSetYieldFunction,  ( WSLINK wslp, WSYieldFunctionObject yf));
#else
WSDECL( wsapi_result,  WSSetMessageHandler, ( WSLINK wslp, WSMessageHandlerObject h));
WSDECL( wsapi_result,  WSSetYieldFunction,  ( WSLINK wslp, WSYieldFunctionObject yf));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE > 1 && WSINTERFACE < 4
#if WSINTERFACE == 3
WSDECL( int,  WSSetYieldFunction0,  ( WSLINK wslp, WSYieldFunctionObject yf, WSLINK cookie));
WSDECL( int,  WSSetMessageHandler0, ( WSLINK wslp, WSMessageHandlerObject func, WSLINK cookie));
#else
WSDECL( wsapi_result,  WSSetYieldFunction0,  ( WSLINK wslp, WSYieldFunctionObject yf, WSLINK cookie));
WSDECL( wsapi_result,  WSSetMessageHandler0, ( WSLINK wslp, WSMessageHandlerObject func, WSLINK cookie));
#endif /* WSINTERFACE == 3 */
#endif /* WSINTERFACE > 1 && WSINTERFACE < 4 */

#if WSINTERFACE >= 3
WSDECL( int, WSDeviceInformation, ( WSLINK wslp, devinfo_selector selector, void* buf, long *buflen));
#else
WSDECL( wsapi_result, WSDeviceInformation, ( WSLINK wslp, devinfo_selector selector, WSPointer buf, longp_st buflen));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL( int,         WSLowLevelDeviceName, (WSLINK wslp, const char **name));
WSDECL( void, WSReleaseLowLevelDeviceName, (WSLINK wslp, const char *name));
#endif


WS_END_EXTERN_C



/*************** Textual interface ***************/


#ifndef _WSGET_H
#define _WSGET_H




#endif /* _WSGET_H */

/* explicitly not protected by _WSGET_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( int,   WSGetNext,          ( WSLINK wslp));
WSDECL( int,   WSGetNextRaw,       ( WSLINK wslp));
WSDECL( int,   WSGetType,          ( WSLINK wslp));
WSDECL( int,   WSGetRawType,       ( WSLINK wslp));
WSDECL( int,   WSGetRawData,       ( WSLINK wslp, unsigned char *data, int size, int *gotp));
WSDECL( int,   WSGetData,          ( WSLINK wslp, char *data, int size, int *gotp));
WSDECL( int,   WSGetArgCount,      ( WSLINK wslp, int *countp));
WSDECL( int,   WSGetRawArgCount,   ( WSLINK wslp, int *countp));
WSDECL( int,   WSBytesToGet,       ( WSLINK wslp, int *leftp));
WSDECL( int,   WSRawBytesToGet,    ( WSLINK wslp, int *leftp));
WSDECL( int,   WSExpressionsToGet, ( WSLINK wslp, int *countp));
#else
WSDECL( wsapi_token,    WSGetNext,          ( WSLINK wslp));
WSDECL( wsapi_token,    WSGetNextRaw,       ( WSLINK wslp));
WSDECL( wsapi_token,    WSGetType,          ( WSLINK wslp));
WSDECL( wsapi_token,    WSGetRawType,       ( WSLINK wslp));
WSDECL( wsapi_result,   WSGetRawData,       ( WSLINK wslp, ucharp_ct data, long_st size, longp_st gotp));
WSDECL( wsapi_result,   WSGetData,          ( WSLINK wslp, charp_ct data, long_st size, longp_st gotp));
WSDECL( wsapi_result,   WSGetArgCount,      ( WSLINK wslp, longp_st countp));
WSDECL( wsapi_result,   WSGetRawArgCount,   ( WSLINK wslp, longp_st countp));
WSDECL( wsapi_result,   WSBytesToGet,       ( WSLINK wslp, longp_st leftp));
WSDECL( wsapi_result,   WSRawBytesToGet,    ( WSLINK wslp, longp_st leftp));
WSDECL( wsapi_result,   WSExpressionsToGet, ( WSLINK wslp, longp_st countp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,   WSNewPacket,        ( WSLINK wslp));
#else
WSDECL( wsapi_result,   WSNewPacket,        ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,   WSTakeLast,         ( WSLINK wslp, int eleft));
#else
WSDECL( wsapi_result,   WSTakeLast,         ( WSLINK wslp, long_st eleft));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,   WSFill,             ( WSLINK wslp));
#else
WSDECL( wsapi_result,   WSFill,             ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C




#ifndef _WSPUT_H
#define _WSPUT_H






#define WSPutExpression is obsolete, use WSPutComposite

#endif /* _WSPUT_H */

/* explicitly not protected by _WSPUT_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( int,   WSPutNext,      ( WSLINK wslp, int tok));
WSDECL( int,   WSPutType,      ( WSLINK wslp, int tok));
WSDECL( int,   WSPutRawSize,   ( WSLINK wslp, int size));
WSDECL( int,   WSPutRawData,   ( WSLINK wslp, const unsigned char *data, int len));
WSDECL( int,   WSPutArgCount,  ( WSLINK wslp, int argc));
WSDECL( int,   WSPutComposite, ( WSLINK wslp, int argc));
WSDECL( int,   WSBytesToPut,   ( WSLINK wslp, int *leftp));
WSDECL( int,   WSEndPacket,    ( WSLINK wslp));
WSDECL( int,   WSFlush,        ( WSLINK wslp));
#else
WSDECL( wsapi_result,   WSPutNext,      ( WSLINK wslp, wsapi_token tok));
WSDECL( wsapi_result,   WSPutType,      ( WSLINK wslp, wsapi__token tok));
WSDECL( wsapi_result,   WSPutRawSize,   ( WSLINK wslp, long_st size));
WSDECL( wsapi_result,   WSPutRawData,   ( WSLINK wslp, kucharp_ct data, long_st len));
WSDECL( wsapi_result,   WSPutArgCount,  ( WSLINK wslp, long_st argc));
WSDECL( wsapi_result,   WSPutComposite, ( WSLINK wslp, long_st argc));
WSDECL( wsapi_result,   WSBytesToPut,   ( WSLINK wslp, longp_st leftp));
WSDECL( wsapi_result,   WSEndPacket,    ( WSLINK wslp));
WSDECL( wsapi_result,   WSFlush,        ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C




#ifndef _WSTK_H
#define _WSTK_H


#define	WSTKOLDINT     'I'		/* 73 Ox49 01001001 */ /* integer leaf node */
#define	WSTKOLDREAL    'R'		/* 82 Ox52 01010010 */ /* real leaf node */


#define	WSTKFUNC    'F'		/* 70 Ox46 01000110 */ /* non-leaf node */

#define	WSTKERROR   (0)		/* bad token */
#define	WSTKERR     (0)		/* bad token */

/* text token bit patterns: 0010x01x --exactly 2 bits worth chosen to make things somewhat readable */
#define WSTK__IS_TEXT( tok) ( (tok & 0x00F6) == 0x0022)


#define	WSTKSTR     '"'         /* 34 0x22 00100010 */
#define	WSTKSYM     '\043'      /* 35 0x23 # 00100011 */ /* octal here as hash requires a trigraph */

#if WSINTERFACE >= 4
#define WSTKOPTSYM  'O'       /* 79 00101010 */
#define WSTKOPTSTR  'Q'       /* 81 01010001 */
#endif

#define	WSTKREAL    '*'         /* 42 0x2A 00101010 */
#define	WSTKINT     '+'         /* 43 0x2B 00101011 */



/* The following defines are for internal use only */
#define	WSTKPCTEND  ']'     /* at end of top level expression */
#define	WSTKAPCTEND '\n'    /* at end of top level expression */
#define	WSTKEND     '\n'
#define	WSTKAEND    '\r'
#define	WSTKSEND    ','

#define	WSTKCONT    '\\'
#define	WSTKELEN    ' '

#define	WSTKNULL    '.'
#define	WSTKOLDSYM  'Y'     /* 89 0x59 01011001 */
#define	WSTKOLDSTR  'S'     /* 83 0x53 01010011 */


typedef unsigned long decoder_mask;
#define	WSTKPACKED	'P'     /* 80 0x50 01010000 */
#define	WSTKARRAY	'A'     /* 65 0x41 01000001 */
#define	WSTKDIM		'D'     /* 68 0x44 01000100 */

#define WSLENGTH_DECODER        ((decoder_mask) 1<<16)
#define WSTKPACKED_DECODER      ((decoder_mask) 1<<17)
#define WSTKARRAY_DECODER	    ((decoder_mask) 1<<18)
#define WSTKMODERNCHARS_DECODER ((decoder_mask) 1<<19)
#if 0
#define WSTKNULLSEQUENCE_DECODER ((decoder_mask) 1<<20)
#else
#define WSTKNULLSEQUENCE_DECODER ((decoder_mask) 0)
#endif
#define WSTKALL_DECODERS (WSLENGTH_DECODER | WSTKPACKED_DECODER | WSTKARRAY_DECODER | WSTKMODERNCHARS_DECODER | WSTKNULLSEQUENCE_DECODER)

#define WSTK_FIRSTUSER '\x30' /* user token */
#define WSTK_LASTUSER  '\x3F'



#endif /* _WSTK_H */



/*************** wsint64 interface ***************/


#ifndef _WSINTEGER64_H
#define _WSINTEGER64_H



#if WSINTERFACE == 3

WS_EXTERN_C

WSDECL(wsint64,              WSInteger64_FromInteger16,    (short a));
WSDECL(wsint64,              WSInteger64_FromInteger32,    (int a));
WSDECL(short,                  WSInteger64_ToInteger16,    (wsint64 a));
WSDECL(int,                    WSInteger64_ToInteger32,    (wsint64 a));

WSDECL(int,                          WSInteger64_Equal,    (wsint64 a, wsint64 b));
WSDECL(int,                       WSInteger64_NotEqual,    (wsint64 a, wsint64 b));
WSDECL(int,                      WSInteger64_LessEqual,    (wsint64 a, wsint64 b));
WSDECL(int,                           WSInteger64_Less,    (wsint64 a, wsint64 b));
WSDECL(int,                   WSInteger64_GreaterEqual,    (wsint64 a, wsint64 b));
WSDECL(int,                        WSInteger64_Greater,    (wsint64 a, wsint64 b));

WSDECL(int,    WSInteger64_GreaterThanLargestInteger32,    (wsint64 a));

WSDECL(wsint64,                     WSInteger64_BitAnd,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                      WSInteger64_BitOr,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                     WSInteger64_BitXor,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                     WSInteger64_BitNot,    (wsint64 a));
WSDECL(wsint64,                        WSInteger64_Neg,    (wsint64 a));

WSDECL(wsint64,                       WSInteger64_Plus,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                  WSInteger64_Increment,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                WSInteger64_Increment32,    (wsint64 a, int b));
WSDECL(wsint64,                   WSInteger64_Subtract,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                  WSInteger64_Decrement,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                WSInteger64_Decrement32,    (wsint64 a, int b));
WSDECL(wsint64,                      WSInteger64_Times,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                     WSInteger64_Divide,    (wsint64 a, wsint64 b));
WSDECL(wsint64,                        WSInteger64_Mod,    (wsint64 a, wsint64 b));

WS_END_EXTERN_C

#endif /* WSINTERFACE == 3 */

#endif /* _WSINTEGER64_H */



/*************** Native C types interface ***************/


#ifndef WSGETNUMBERS_HPP
#define WSGETNUMBERS_HPP







#if WSINTERFACE < 3
#define WSGetReal WSGetDouble
#endif

#endif /* WSGETNUMBERS_HPP */


/* explicitly not protected by WSGETNUMBERS_HPP in case WSDECL is redefined for multiple inclusion */


WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( int,   WSGetBinaryNumber,  ( WSLINK wslp, void *np, long type));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetShortInteger   - WSGetInteger16
WSGetInteger        - WSGetInteger32
WSGetLongInteger    - WSGetInteger64 for 64-bit integers or WSGetInteger32 for 32-bit integers
*/
WSDECL( int,   WSGetShortInteger,  ( WSLINK wslp, short *hp));
WSDECL( int,   WSGetInteger,       ( WSLINK wslp, int *ip));
WSDECL( int,   WSGetLongInteger,   ( WSLINK wslp, long *lp));


WSDECL( int,   WSGetInteger16,  ( WSLINK wslp, short *hp));
WSDECL( int,   WSGetInteger32,  ( WSLINK wslp, int *ip));
WSDECL( int,   WSGetInteger64,  ( WSLINK wslp, wsint64 *wp));
#else
WSDECL( wsapi_result,   WSGetBinaryNumber,  ( WSLINK wslp, voidp_ct np, long type));
WSDECL( wsapi_result,   WSGetShortInteger,  ( WSLINK wslp, shortp_nt hp));
WSDECL( wsapi_result,   WSGetInteger,       ( WSLINK wslp, intp_nt ip));
WSDECL( wsapi_result,   WSGetLongInteger,   ( WSLINK wslp, longp_nt lp));
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE >= 4
WSDECL(int, WSGetInteger8, (WSLINK wslp, unsigned char *cp));
#endif


#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetFloat      - WSGetReal32
WSGetDouble     - WSGetReal64
WSGetReal       - WSGetReal64
WSGetLongDouble - WSGetReal128
*/
WSDECL( int,   WSGetFloat,         ( WSLINK wslp, float *fp));
WSDECL( int,   WSGetDouble,        ( WSLINK wslp, double *dp));
WSDECL( int,   WSGetReal,          ( WSLINK wslp, double *dp));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSGetLongDouble,    ( WSLINK wslp, wsextended_double *xp));
#endif


WSDECL( int,   WSGetReal32,         ( WSLINK wslp, float *fp));
WSDECL( int,   WSGetReal64,        ( WSLINK wslp, double *dp));
WSDECL( int,   WSGetReal128,          ( WSLINK wslp, wsextended_double *dp));
#else
WSDECL( wsapi_result,   WSGetFloat,         ( WSLINK wslp, floatp_nt fp));
WSDECL( wsapi_result,   WSGetDouble,        ( WSLINK wslp, doublep_nt dp));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSGetLongDouble,    ( WSLINK wslp, extendedp_nt xp));
#endif
#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C





#ifndef WSGETSTRINGS_HPP
#define WSGETSTRINGS_HPP







#endif /* WSGETSTRINGS_HPP */

/* explicitly not protected by WSGETSTRINGS_HPP in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 WSGet16BitCharacters has been deprecated.  Use the suggested function in its
place:

WSGet16BitCharacters   - WSGetUCS2Characters
*/

#if WSINTERFACE == 3
WSDECL( int,   WSGet16BitCharacters,  ( WSLINK wslp, long *chars_left, unsigned short *buf, long cardof_buf, long *got));
#endif

WSDECL( int,   WSGet8BitCharacters,   ( WSLINK wslp, long *chars_left, unsigned char *buf, long cardof_buf, long *got, long missing));
WSDECL( int,   WSGet7BitCharacters,   ( WSLINK wslp, long *chars_left, char *buf, long cardof_buf, long *got));

WSDECL( int,   WSGetUCS2Characters,   ( WSLINK wslp, int *chars_left, unsigned short *buf, int cardof_buf, int *got));
WSDECL( int,   WSGetUTF8Characters,   ( WSLINK wslp, int *chars_left, unsigned char *buf, int cardof_buf, int *got));
WSDECL( int,   WSGetUTF16Characters,  ( WSLINK wslp, int *chars_left, unsigned short *buf, int cardof_buf, int *got));
WSDECL( int,   WSGetUTF32Characters,  ( WSLINK wslp, int *chars_left, unsigned int *buf, int cardof_buf, int *got));

/*
As of WSINTERFACE 3 WSGetUnicodeString has been deprecated.  Use the suggested function in its
place:

WSGetUnicodeString - WSGetUCS2String
*/

#if WSINTERFACE == 3
WSDECL( int,   WSGetUnicodeString,    ( WSLINK wslp, const unsigned short **sp, long *lenp));
#endif

WSDECL( int,   WSGetByteString,       ( WSLINK wslp, const unsigned char **sp, int *lenp, long missing));
WSDECL( int,   WSGetString,           ( WSLINK wslp, const char **sp));

WSDECL( int,   WSGetUCS2String,       ( WSLINK wslp, const unsigned short **sp, int *lenp));
WSDECL( int,   WSGetUTF8String,       ( WSLINK wslp, const unsigned char **sp, int *bytes, int *chars));
WSDECL( int,   WSGetUTF16String,      ( WSLINK wslp, const unsigned short **sp, int *ncodes, int *chars));
WSDECL( int,   WSGetUTF32String,      ( WSLINK wslp, const unsigned int **sp, int *len));
#else
WSDECL( wsapi_result,   WSGet16BitCharacters,  ( WSLINK wslp, longp_st chars_left, ushortp_ct buf, long_st cardof_buf, longp_st got));
WSDECL( wsapi_result,   WSGet8BitCharacters,   ( WSLINK wslp, longp_st chars_left, ucharp_ct  buf, long_st cardof_buf, longp_st got, long missing));
WSDECL( wsapi_result,   WSGet7BitCharacters,   ( WSLINK wslp, longp_st chars_left, charp_ct   buf, long_st cardof_buf, longp_st got));

WSDECL( wsapi_result,   WSGetUnicodeString,    ( WSLINK wslp, kushortpp_ct sp, longp_st lenp));
WSDECL( wsapi_result,   WSGetByteString,       ( WSLINK wslp, kucharpp_ct  sp, longp_st lenp, long missing));
WSDECL( wsapi_result,   WSGetString,           ( WSLINK wslp, kcharpp_ct   sp));
#endif /* WSINTERFACE >= 3 */

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 WSGetUnicodeString0 has been deprecated.  Use the suggested function in its
place:

WSGetUnicodeString0 - WSGetUCS2String0
*/

#if WSINTERFACE >= 4

WSDECL( int,   WSGetNumberAsByteString,      ( WSLINK wslp, const unsigned char **sp, long *lenp, long missing));
WSDECL( int,   WSGetNumberAsString,          ( WSLINK wslp, const char **sp));

WSDECL( int,   WSGetNumberAsUCS2String,      ( WSLINK wslp, const unsigned short **sp, int *lenp));
WSDECL( int,   WSGetNumberAsUTF8String,      ( WSLINK wslp, const unsigned char **sp, int *bytes, int *chars));
WSDECL( int,   WSGetNumberAsUTF16String,     ( WSLINK wslp, const unsigned short **sp, int *ncodes, int *chars));
WSDECL( int,   WSGetNumberAsUTF32String,     ( WSLINK wslp, const unsigned int **sp, int *lenp));
#else
WSDECL( int,   WSGetUnicodeString0,   ( WSLINK wslp, const unsigned short **sp, long *lenp));
WSDECL( int,   WSGetByteString0,      ( WSLINK wslp, const unsigned char **sp, long *lenp, long missing));
WSDECL( int,   WSGetString0,          ( WSLINK wslp, const char **sp));

WSDECL( int,   WSGetUCS2String0,      ( WSLINK wslp, const unsigned short **sp, int *lenp));
WSDECL( int,   WSGetUTF8String0,      ( WSLINK wslp, const unsigned char **sp, int *bytes, int *chars));
WSDECL( int,   WSGetUTF16String0,     ( WSLINK wslp, const unsigned short **sp, int *ncodes, int *chars));
WSDECL( int,   WSGetUTF32String0,     ( WSLINK wslp, const unsigned int **sp, int *lenp));
#endif /* WSINTERFACE >= 4 */
#else
#if WSINTERFACE > 1
WSDECL( wsapi_result,   WSGetUnicodeString0,   ( WSLINK wslp, kushortpp_ct sp, longp_st lenp));
WSDECL( wsapi_result,   WSGetByteString0,      ( WSLINK wslp, kucharpp_ct  sp, longp_st lenp, long missing));
WSDECL( wsapi_result,   WSGetString0,          ( WSLINK wslp, kcharpp_ct   sp));
#endif
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSDisownUnicodeString - WSReleaseUCS2String 
WSDisownByteString    - WSReleaseByteString
WSDisownString        - WSReleaseString
*/

#if WSINTERFACE == 3
WSDECL( void,           WSDisownUnicodeString, ( WSLINK wslp, const unsigned short *s,   long len));
WSDECL( void,           WSDisownByteString,    ( WSLINK wslp, const unsigned char * s,   long len));
WSDECL( void,           WSDisownString,        ( WSLINK wslp, const char *s));
#endif

WSDECL( void,           WSReleaseUCS2String,   ( WSLINK wslp, const unsigned short *s,   int len));
WSDECL( void,           WSReleaseUTF8String,   ( WSLINK wslp, const unsigned char *s, int len));
WSDECL( void,           WSReleaseUTF16String,  ( WSLINK wslp, const unsigned short *s, int len));
WSDECL( void,           WSReleaseUTF32String,  ( WSLINK wslp, const unsigned int *s, int len));
WSDECL( void,           WSReleaseByteString,   ( WSLINK wslp, const unsigned char * s,   int len));
WSDECL( void,           WSReleaseString,       ( WSLINK wslp, const char *s));

#if WSINTERFACE <= 3
WSDECL( int,   WSCheckString,   ( WSLINK wslp, const char *name));
#endif

#if WSINTERFACE >= 4
WSDECL( int,    WSTestString,      ( WSLINK wslp, const char *name));
WSDECL( int,    WSTestUCS2String,  ( WSLINK wslp, const unsigned short *name, int length));
WSDECL( int,    WSTestUTF8String,  ( WSLINK wslp, const unsigned char *name, int length));
WSDECL( int,    WSTestUTF16String, ( WSLINK wslp, const unsigned short *name, int length));
WSDECL( int,    WSTestUTF32String, ( WSLINK wslp, const unsigned int *name, int length));
#endif // WSINTERFACE >= 4

#else


WSDECL( void,           WSDisownUnicodeString, ( WSLINK wslp, kushortp_ct s,   long_st len));
WSDECL( void,           WSDisownByteString,    ( WSLINK wslp, kucharp_ct  s,   long_st len));
WSDECL( void,           WSDisownString,        ( WSLINK wslp, kcharp_ct   s));

WSDECL( wsapi_result,   WSCheckString,   ( WSLINK wslp, kcharp_ct name));

#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C





#ifndef WSGETSYMBOLS_HPP
#define WSGETSYMBOLS_HPP







#endif /* WSGETSYMBOLS_HPP */

/* explicitly not protected by WSGETSYMBOLS_HPP in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3

/*
As of WSINTERFACE 3 WSGetUnicodeSymbol has been deprecated.  Use the suggested function in its
place:

WSGetUnicodeSymbol - WSGetUCS2Symbol
*/

#if WSINTERFACE == 3
WSDECL( int,   WSGetUnicodeSymbol,    ( WSLINK wslp, const unsigned short **sp, long *lenp));
#endif

WSDECL( int,   WSGetByteSymbol,       ( WSLINK wslp, const unsigned char ** sp, int *lenp, long missing));
WSDECL( int,   WSGetSymbol,           ( WSLINK wslp, const char **          sp));

WSDECL( int,   WSGetUCS2Symbol,       ( WSLINK wslp, const unsigned short **sp, int *lenp));
WSDECL( int,   WSGetUTF8Symbol,       ( WSLINK wslp, const unsigned char **sp, int *bytes, int *chars));
WSDECL( int,   WSGetUTF16Symbol,      ( WSLINK wslp, const unsigned short **sp, int *ncodes, int *chars));
WSDECL( int,   WSGetUTF32Symbol,      ( WSLINK wslp, const unsigned int **sp, int *lenp));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:
WSDisownUnicodeSymbol - WSReleaseUCS2Symbol
WSDisownByteSymbol    - WSReleaseByteSymbol
WSDisownSymbol        - WSReleaseSymbol
*/

#if WSINTERFACE == 3
WSDECL( void,           WSDisownUnicodeSymbol, ( WSLINK wslp, const unsigned short *s,   long len));
WSDECL( void,           WSDisownByteSymbol,    ( WSLINK wslp, const unsigned char * s,   long len));
WSDECL( void,           WSDisownSymbol,        ( WSLINK wslp, const char *s));
#endif

WSDECL( void,           WSReleaseUCS2Symbol,   ( WSLINK wslp, const unsigned short *s,   int len));
WSDECL( void,           WSReleaseUTF8Symbol,   ( WSLINK wslp, const unsigned char *s, int len));
WSDECL( void,           WSReleaseUTF16Symbol,  ( WSLINK wslp, const unsigned short *s, int len));
WSDECL( void,           WSReleaseUTF32Symbol,  ( WSLINK wslp, const unsigned int *s, int len));
WSDECL( void,           WSReleaseByteSymbol,   ( WSLINK wslp, const unsigned char * s,   int len));
WSDECL( void,           WSReleaseSymbol,       ( WSLINK wslp, const char *s));

#if WSINTERFACE == 3
WSDECL( int,            WSCheckSymbol,         ( WSLINK wslp, const char *name));
#endif

#if WSINTERFACE >= 4
WSDECL( int,            WSTestSymbol,          ( WSLINK wslp, const char *name));
WSDECL( int,            WSTestUCS2Symbol,      ( WSLINK wslp, const unsigned short *name, int length));
WSDECL( int,            WSTestUTF8Symbol,      ( WSLINK wslp, const unsigned char *name, int length));
WSDECL( int,            WSTestUTF16Symbol,     ( WSLINK wslp, const unsigned short *name, int length));
WSDECL( int,            WSTestUTF32Symbol,     ( WSLINK wslp, const unsigned int *name, int length));
#endif

WSDECL( int,            WSGetFunction,         ( WSLINK wslp, const char **sp, int *countp));

#if WSINTERFACE >= 4
WSDECL( int,            WSGetUCS2Function,     ( WSLINK wslp, const unsigned short **sp, int *length, int *countp));
WSDECL( int,            WSGetUTF8Function,     ( WSLINK wslp, const unsigned char **sp, int *length, int *countp));
WSDECL( int,            WSGetUTF16Function,    ( WSLINK wslp, const unsigned short **sp, int *length, int *countp));
WSDECL( int,            WSGetUTF32Function,    ( WSLINK wslp, const unsigned int **sp, int *length, int *countp));
#endif

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSCheckFunction             - WSTestHead
WSCheckFunctionWithArgCount - WSTestHead
*/

WSDECL( int,   WSCheckFunction, ( WSLINK wslp, const char *s, long *countp));
WSDECL( int,   WSCheckFunctionWithArgCount, ( WSLINK wslp, const char *s, long *countp));

WSDECL( int,   WSTestHead,      ( WSLINK wslp, const char *s, int *countp));

#if WSINTERFACE >= 4

/*
For a limited time convenience define the following:
*/

WSDECL( int,      WSTestHeadWithArgCount, (WSLINK wslp, const char *s, int *countp));
WSDECL( int,  WSTestUCS2HeadWithArgCount, (WSLINK wslp, const unsigned short *s, int length, int *countp));
WSDECL( int, WSTestUTF16HeadWithArgCount, (WSLINK wslp, const unsigned short *s, int length, int *countp));
WSDECL( int, WSTestUTF32HeadWithArgCount, (WSLINK wslp, const unsigned int *s, int length, int *countp));
WSDECL( int,  WSTestUTF8HeadWithArgCount, (WSLINK wslp, const unsigned char *s, int length, int *countp));

WSDECL( int,   WSTestUCS2Head,  ( WSLINK wslp, const unsigned short *s, int length, int *countp));
WSDECL( int,   WSTestUTF8Head,  ( WSLINK wslp, const unsigned char *s, int length, int *countp));
WSDECL( int,   WSTestUTF16Head, ( WSLINK wslp, const unsigned short *s, int length, int *countp));
WSDECL( int,   WSTestUTF32Head, ( WSLINK wslp, const unsigned int *s, int length, int *countp));
#endif // WSINTERFACE >= 3

#else
WSDECL( wsapi_result,   WSGetUnicodeSymbol,    ( WSLINK wslp, kushortpp_ct sp, longp_st lenp));
WSDECL( wsapi_result,   WSGetByteSymbol,       ( WSLINK wslp, kucharpp_ct  sp, longp_st lenp, long missing));
WSDECL( wsapi_result,   WSGetSymbol,           ( WSLINK wslp, kcharpp_ct   sp));

WSDECL( void,           WSDisownUnicodeSymbol, ( WSLINK wslp, kushortp_ct s,   long_st len));
WSDECL( void,           WSDisownByteSymbol,    ( WSLINK wslp, kucharp_ct  s,   long_st len));
WSDECL( void,           WSDisownSymbol,        ( WSLINK wslp, kcharp_ct   s));

WSDECL( wsapi_result,   WSCheckSymbol,   ( WSLINK wslp, kcharp_ct name));
WSDECL( wsapi_result,   WSGetFunction,   ( WSLINK wslp, kcharpp_ct sp, longp_st countp));
WSDECL( wsapi_result,   WSCheckFunction, ( WSLINK wslp, kcharp_ct s, longp_st countp));
WSDECL( wsapi_result,   WSCheckFunctionWithArgCount, ( WSLINK wslp, kcharp_ct s, longp_st countp));

#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C







#ifndef WSPUTNUMBERS_HPP
#define WSPUTNUMBERS_HPP







#if WSINTERFACE < 3
#define WSPutReal WSPutDouble
#endif

#endif /* WSPUTNUMBERS_HPP */


/* explicitly not protected by WSPUTNUMBERS_HPP in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( int,   WSPutBinaryNumber,  ( WSLINK wslp, void *np, long type));

/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutShortInteger  - WSPutInteger16
WSPutInteger       - WSPutInteger32
WSPutLongInteger   - WSPutInteger64 for 64-bit integers or WSPutInteger32 for 32-bit integers.
*/
WSDECL( int,   WSPutShortInteger,  ( WSLINK wslp, int h));
WSDECL( int,   WSPutInteger,       ( WSLINK wslp, int i));
WSDECL( int,   WSPutLongInteger,   ( WSLINK wslp, long l));

WSDECL( int,   WSPutInteger16,     ( WSLINK wslp, int h));
WSDECL( int,   WSPutInteger32,     ( WSLINK wslp, int i));
WSDECL( int,   WSPutInteger64,     ( WSLINK wslp, wsint64 w));
#else
WSDECL( wsapi_result,   WSPutBinaryNumber,  ( WSLINK wslp, voidp_ct np, long type));
WSDECL( wsapi_result,   WSPutShortInteger,  ( WSLINK wslp, int_nt h));
WSDECL( wsapi_result,   WSPutInteger,       ( WSLINK wslp, int_nt i));
WSDECL( wsapi_result,   WSPutLongInteger,   ( WSLINK wslp, long_nt l));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL( int, WSPutInteger8,   (WSLINK wslp, unsigned char i));
#endif

#if WSINTERFACE >= 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutFloat      - WSPutReal32
WSPutDouble     - WSPutReal64
WSPutReal       - WSPutReal64
WSPutLongDouble - WSPutReal128
*/
WSDECL( int,   WSPutFloat,         ( WSLINK wslp, double f));
WSDECL( int,   WSPutDouble,        ( WSLINK wslp, double d));
WSDECL( int,   WSPutReal,          ( WSLINK wslp, double d));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSPutLongDouble,    ( WSLINK wslp, wsextended_double x));
#endif

WSDECL( int,   WSPutReal32,         ( WSLINK wslp, double f));
WSDECL( int,   WSPutReal64,         ( WSLINK wslp, double d));
WSDECL( int,   WSPutReal128,        ( WSLINK wslp, wsextended_double x));
#else
WSDECL( wsapi_result,   WSPutFloat,         ( WSLINK wslp, double_nt f));
WSDECL( wsapi_result,   WSPutDouble,        ( WSLINK wslp, double_nt d));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSPutLongDouble,   ( WSLINK wslp, extended_nt x));
#endif
#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C





#ifndef WSPUTSTRINGS_HPP
#define WSPUTSTRINGS_HPP








#endif /* WSPUTSTRINGS_HPP */

/* explicitly not protected by WSPUTSTRINGS_HPP in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 WSPut16BitCharacters has been deprecated.  Use the suggested function in its
place:

WSPut16BitCharacters   - WSPutUCS2Characters
*/

#if WSINTERFACE == 3
WSDECL( int,   WSPut16BitCharacters, ( WSLINK wslp, long chars_left, const unsigned short *codes, long ncodes));
#endif

WSDECL( int,   WSPut8BitCharacters,  ( WSLINK wslp, long chars_left, const unsigned char *bytes, long nbytes));
WSDECL( int,   WSPut7BitCount,       ( WSLINK wslp, long count, long size));
WSDECL( int,   WSPut7BitCharacters,  ( WSLINK wslp, long chars_left, const char *bytes, long nbytes, long nchars_now));

WSDECL( int,   WSPutUCS2Characters,  ( WSLINK wslp, int chars_left, const unsigned short *codes, int ncodes));
WSDECL( int,   WSPutUTF8Characters,  ( WSLINK wslp, int chars_left, const unsigned char *codes, int ncodes));
WSDECL( int,   WSPutUTF16Characters, ( WSLINK wslp, int chars_left, const unsigned short *codes, int ncodes));
WSDECL( int,   WSPutUTF32Characters, ( WSLINK wslp, int chars_left, const unsigned int *codes, int ncodes));

/*
As of WSINTERFACE 3 WSPutUnicodeString has been deprecated.  Use the suggested function in its
place:

WSPutUnicodeString - WSPutUCS2String
*/

#if WSINTERFACE == 3
WSDECL( int,   WSPutUnicodeString, ( WSLINK wslp, const unsigned short *s, long len));
#endif

WSDECL( int,   WSPutByteString,    ( WSLINK wslp, const unsigned char *s, long len));
WSDECL( int,   WSPutString,        ( WSLINK wslp, const char *s));

WSDECL( int,   WSPutUCS2String,    ( WSLINK wslp, const unsigned short *s, int len));
WSDECL( int,   WSPutUTF8String,    ( WSLINK wslp, const unsigned char *s, int len));
WSDECL( int,   WSPutUTF16String,   ( WSLINK wslp, const unsigned short *s, int len));
WSDECL( int,   WSPutUTF32String,   ( WSLINK wslp, const unsigned int *s, int len));
#else
WSDECL( wsapi_result,   WSPut16BitCharacters, ( WSLINK wslp, long_st chars_left, kushortp_ct codes, long_st ncodes));
WSDECL( wsapi_result,   WSPut8BitCharacters,  ( WSLINK wslp, long_st chars_left, kucharp_ct bytes, long_st nbytes));
WSDECL( wsapi_result,   WSPut7BitCount,       ( WSLINK wslp, long_st count, long_st size));
WSDECL( wsapi_result,   WSPut7BitCharacters,  ( WSLINK wslp, long_st chars_left, kcharp_ct bytes, long_st nbytes, long_st nchars_now));

WSDECL( wsapi_result,   WSPutUnicodeString, ( WSLINK wslp, kushortp_ct s, long_st len));
WSDECL( wsapi_result,   WSPutByteString,    ( WSLINK wslp, kucharp_ct  s, long_st len));
WSDECL( wsapi_result,   WSPutString,        ( WSLINK wslp, kcharp_ct   s));
#endif /* WSINTERFACE >= 3 */

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 WSPutRealUnicodeString0 has been deprecated.  Use the suggested function in its
place:

WSPutRealUnicodeString0 - WSPutRealUCS2String0
*/

#if WSINTERFACE >= 4
WSDECL( int,   WSPutRealNumberAsString,        ( WSLINK wslp, const char *s));
WSDECL( int,   WSPutRealNumberAsByteString,    ( WSLINK wslp, const unsigned char *s));
WSDECL( int,   WSPutRealNumberAsUCS2String,    ( WSLINK wslp, const unsigned short *s));
WSDECL( int,   WSPutRealNumberAsUTF8String,    ( WSLINK wslp, const unsigned char *s, int nbytes));
WSDECL( int,   WSPutRealNumberAsUTF16String,   ( WSLINK wslp, const unsigned short *s, int ncodes));
WSDECL( int,   WSPutRealNumberAsUTF32String,   ( WSLINK wslp, const unsigned int *s, int nchars));
#endif /* WSINTERFACE >= 4 */

#if WSINTERFACE == 3
WSDECL( int,   WSPutRealUnicodeString0, ( WSLINK wslp, unsigned short *s));
WSDECL( int,   WSPutRealByteString0,    ( WSLINK wslp, unsigned char *s));

WSDECL( int,   WSPutRealUCS2String0,    ( WSLINK wslp, unsigned short *s));
WSDECL( int,   WSPutRealUTF8String0,    ( WSLINK wslp, unsigned char *s, int nbytes));
WSDECL( int,   WSPutRealUTF16String0,   ( WSLINK wslp, unsigned short *s, int ncodes));
WSDECL( int,   WSPutRealUTF32String0,   ( WSLINK wslp, unsigned int *s, int nchars));
#endif /* WSINTERFACE == 3 */
#else
#if WSINTERFACE > 1
WSDECL( wsapi_result,   WSPutRealUnicodeString0, ( WSLINK wslp, ushortp_ct s));
WSDECL( wsapi_result,   WSPutRealByteString0,    ( WSLINK wslp, ucharp_ct  s));
#endif
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,   WSPutSize,          ( WSLINK wslp, int size));
WSDECL( int,   WSPutData,          ( WSLINK wslp, const char *buff, int len));
#else
WSDECL( wsapi_result,   WSPutSize, ( WSLINK wslp, long_st size));
WSDECL( wsapi_result,   WSPutData, ( WSLINK wslp, kcharp_ct buff, long_st len));
#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C





#ifndef WSPUTSYMBOLS_HPP
#define WSPUTSYMBOLS_HPP








#endif /* WSPUTSYMBOLS_HPP */

/* explicitly not protected by WSPUTSYMBOLS_HPP in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 WSPutUnicodeSymbol has been deprecated.  Use the suggested function in its
place:

WSPutUnicodeSymbol - WSPutUCS2Symbol
*/

#if WSINTERFACE == 3
WSDECL( int,   WSPutUnicodeSymbol, ( WSLINK wslp, const unsigned short *s, long len));
#endif

WSDECL( int,   WSPutByteSymbol,    ( WSLINK wslp, const unsigned char *s, long len));
WSDECL( int,   WSPutSymbol,        ( WSLINK wslp, const char *s));

WSDECL( int,   WSPutUCS2Symbol,    ( WSLINK wslp, const unsigned short *s, int len));

WSDECL( int,   WSPutUTF8Symbol,    ( WSLINK wslp, const unsigned char *s, int len));
WSDECL( int,   WSPutUTF16Symbol,   ( WSLINK wslp, const unsigned short *s, int len));
WSDECL( int,   WSPutUTF32Symbol,   ( WSLINK wslp, const unsigned int *s, int len));


WSDECL( int,   WSPutFunction,      ( WSLINK wslp, const char *s, int argc));

#if WSINTERFACE >= 4
WSDECL( int,   WSPutUCS2Function,  ( WSLINK wslp, const unsigned short *s, int length, int argn));
WSDECL( int,   WSPutUTF8Function,  ( WSLINK wslp, const unsigned char *s, int length, int argn));
WSDECL( int,   WSPutUTF16Function, ( WSLINK wslp, const unsigned short *s, int length, int argn));
WSDECL( int,   WSPutUTF32Function, ( WSLINK wslp, const unsigned int *s, int length, int argn));
#endif // WSINTERFACE >= 4

#else
WSDECL( wsapi_result,   WSPutUnicodeSymbol, ( WSLINK wslp, kushortp_ct s, long_st len));
WSDECL( wsapi_result,   WSPutByteSymbol,    ( WSLINK wslp, kucharp_ct  s, long_st len));
WSDECL( wsapi_result,   WSPutSymbol,        ( WSLINK wslp, kcharp_ct   s));
WSDECL( wsapi_result,   WSPutFunction,      ( WSLINK wslp, kcharp_ct s, long_st argc));
#endif /* WSINTERFACE >= 3 */

WS_END_EXTERN_C





#ifndef _WSSTRING_H
#define _WSSTRING_H










#define MAX_BYTES_PER_OLD_CHARACTER 3
#if WSINTERFACE < 3
#define MAX_BYTES_PER_NEW_CHARACTER 6
#else
#define MAX_BYTES_PER_NEW_CHARACTER 10
#endif

#define WS_MAX_BYTES_PER_CHARACTER MAX_BYTES_PER_NEW_CHARACTER

/* for source code compatibility with earlier versions of WSTP */

typedef struct {
#if WSINTERFACE >= 3
	const char *str;
	const char *end;
#else
	kcharp_ct str;
	kcharp_ct end;
#endif /* WSINTERFACE >= 3 */
} WSStringPosition;

typedef WSStringPosition FAR * WSStringPositionPointer;

#define WSStringFirstPos(s,pos) WSStringFirstPosFun( s, &(pos))

#define WSforString( s, pos) \
	for( WSStringFirstPos(s,pos); WSStringCharacter( (pos).str, (pos).end) >= 0; WSNextCharacter(&(pos).str, (pos).end))

#define WSStringChar( pos) WSStringCharacter( (pos).str, (pos).end)

#define WSPutCharToString WSConvertCharacter


/* for internal use only */

typedef struct {
#if WSINTERFACE >= 3
	unsigned char *cc;
	int  mode;
	int  more;
	unsigned char *head;
#else
	ucharp_ct cc;
	int_ct  mode;
	int_ct  more;
	ucharp_ct head;
#endif /* WSINTERFACE >= 3 */
} WSOldStringPosition;

typedef WSOldStringPosition FAR * WSOldStringPositionPointer;


#define WSOldforString( s, pos) \
  for ( WSOldStringFirstPos( s, pos); (pos).more; WSOldStringNextPos( pos))

#define WSOldStringChar(pos) \
  ( ((pos).mode <= 1) ? (uint_ct)(*(ucharp_ct)((pos).cc)) : WSOldStringCharFun( &pos) )


#define WSOldStringFirstPos(s,pos) WSOldStringFirstPosFun( s, &(pos))

#define WSOldStringNextPos(pos)  ( \
	((pos).mode == 0) \
		? ((*(*(pos).cc ? ++(pos).cc : (pos).cc) ? 0 : ((pos).more = 0)), (pos).cc) \
		: WSOldStringNextPosFun( &pos) )

#endif /* _WSSTRING_H */




WS_EXTERN_C
/* assumes *startp aligned on char boundary, if n == -1, returns ~(char_count) */

#if WSINTERFACE >= 3
WSDECL( long, WSCharacterOffset,           ( const char **startp, const char *end, long n));
WSDECL( long, WSStringCharacter,           ( const char * start,  const char *end));
WSDECL( long, WSNextCharacter,             ( const char **startp, const char *end));
#else
WSDECL( long, WSCharacterOffset,           ( kcharpp_ct startp, kcharp_ct end, long n));
WSDECL( long, WSStringCharacter,           ( kcharp_ct  start,  kcharp_ct end));
WSDECL( long, WSNextCharacter,             ( kcharpp_ct startp, kcharp_ct end));
#endif /* WSINTERFACE >= 3 */


#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif
#if WSINTERFACE > 1
#if WSINTERFACE >= 4
WSDECL( long, WSNextCharacterFromStringWithLength, (const char *str, long *indexp, long len));
#else
#if WSINTERFACE == 3
WSDECL( long, WSNextCharacter0,            ( const char *str, long *indexp, long len));
#else
WSDECL( long, WSNextCharacter0,            ( kcharp_ct str, longp_ct indexp, long len));
#endif /* WSINTERFACE == 3 */
#endif // WSINTERFACE >= 4
#endif // WSINTERFACE > 1

#if WSINTERFACE >= 3
WSDECL( long, WSConvertNewLine,            ( char **sp));
WSDECL( long, WSConvertCharacter,          ( unsigned long ch, char **sp));
WSDECL( long, WSConvertByteString,         ( unsigned char *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertByteStringNL,       ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));
#if WSINTERFACE == 3
WSDECL( long, WSConvertUnicodeString,      ( unsigned short *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertUnicodeStringNL,    ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl));
#endif
WSDECL( long, WSConvertDoubleByteString,   ( unsigned char *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertDoubleByteStringNL, ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));

WSDECL( long, WSConvertUCS2String,         ( unsigned short *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertUCS2StringNL,       ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl));
WSDECL( long, WSConvertUTF8String,         ( unsigned char *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertUTF8StringNL,       ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));
WSDECL( long, WSConvertUTF16String,        ( unsigned short *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertUTF16StringNL,      ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl));
WSDECL( long, WSConvertUTF32String,        ( unsigned int *codes, long len, char **strp, char *str_end));
WSDECL( long, WSConvertUTF32StringNL,      ( unsigned int *codes, long len, char **strp, char *str_end, unsigned long nl));



#else
WSDECL( long, WSConvertNewLine,            ( charpp_ct sp));
WSDECL( long, WSConvertCharacter,          ( ulong_ct ch, charpp_ct sp));
WSDECL( long, WSConvertByteString,         ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
WSDECL( long, WSConvertByteStringNL,       ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));
WSDECL( long, WSConvertUnicodeString,      ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end));
WSDECL( long, WSConvertUnicodeStringNL,    ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));
WSDECL( long, WSConvertDoubleByteString,   ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
WSDECL( long, WSConvertDoubleByteStringNL, ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));
#endif /* WSINTERFACE >= 3 */







/* for source code compatibility with earlier versions of WSTP */
#if WSINTERFACE >= 3
WSDECL( const char *,     WSStringFirstPosFun,  ( const char *s, WSStringPositionPointer p));
#else
WSDECL( kcharp_ct,        WSStringFirstPosFun,  ( kcharp_ct s, WSStringPositionPointer p));
#endif /* WSINTERFACE >= 3 */

/* for internal use only */
#if WSINTERFACE >= 3
WSDECL( int,                WSOldPutCharToString,      ( unsigned int ch, char **sp));
WSDECL( unsigned char *,    WSOldStringNextPosFun,     ( WSOldStringPositionPointer p));
WSDECL( unsigned char *,    WSOldStringFirstPosFun,    ( char *s, WSOldStringPositionPointer p));
WSDECL( unsigned int,       WSOldStringCharFun,        ( WSOldStringPositionPointer p));
WSDECL( long,               WSOldConvertByteString,    ( unsigned char *codes, long len, char **strp, char *str_end));

#if WSINTERFACE == 3
WSDECL( long,               WSOldConvertUnicodeString, ( unsigned short *codes, long len, char **strp, char *str_end));
#endif

WSDECL( long,               WSOldConvertUCS2String,    ( unsigned short *codes, long len, char **strp, char *str_end));
#else
WSDECL( wsapi_result,       WSOldPutCharToString,      ( uint_ct ch, charpp_ct sp));
WSDECL( ucharp_ct,          WSOldStringNextPosFun,     ( WSOldStringPositionPointer p));
WSDECL( ucharp_ct,          WSOldStringFirstPosFun,    ( charp_ct s, WSOldStringPositionPointer p));
WSDECL( uint_ct,            WSOldStringCharFun,        ( WSOldStringPositionPointer p));
WSDECL( long,               WSOldConvertByteString,    ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
WSDECL( long,               WSOldConvertUnicodeString, ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end));
#endif /* WSINTERFACE >= 3 */

/* Internal functions */
WSDECL( long, WSCharOffset,           ( const char **startp, const char *end, long n, int more));
WSDECL( long, WSNextChar,             ( const char **startp, const char *end, int more, int useSurrogates, int *wasSurrogatePair));


WS_END_EXTERN_C




#ifndef _WSCAPUT_H
#define _WSCAPUT_H






#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter FAR * array_meterp;
typedef array_meterp FAR * array_meterpp;
#endif


#if WSINTERFACE < 3
#define WSPutRealArray WSPutDoubleArray
#endif

#endif /* _WSCAPUT_H */


/* explicitly not protected by _WSCAPUT_H in case WSDECL is redefined for multiple inclusion */

/*bugcheck: bugcheck need FAR here */
WS_EXTERN_C
#if WSINTERFACE >= 3
WSDECL( int,   WSPutArray,                  ( WSLINK wslp, array_meterp meterp));
#else
WSDECL( wsapi_result,   WSPutArray,                  ( WSLINK wslp, array_meterp meterp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL( int,   WSPutBinaryNumberArrayData,  ( WSLINK wslp, array_meterp meterp, const void *         datap, long count, long type));
WSDECL( int,   WSPutByteArrayData,          ( WSLINK wslp, array_meterp meterp, const unsigned char *datap, long count));
WSDECL( int,   WSPutShortIntegerArrayData,  ( WSLINK wslp, array_meterp meterp, const short *        datap, long count));
WSDECL( int,   WSPutIntegerArrayData,       ( WSLINK wslp, array_meterp meterp, const int *          datap, long count));
WSDECL( int,   WSPutLongIntegerArrayData,   ( WSLINK wslp, array_meterp meterp, const long *         datap, long count));

WSDECL( int,   WSPutInteger8ArrayData,      ( WSLINK wslp, array_meterp meterp, const unsigned char * datap, int count));
WSDECL( int,   WSPutInteger16ArrayData,     ( WSLINK wslp, array_meterp meterp, const short *        datap, int count));
WSDECL( int,   WSPutInteger32ArrayData,     ( WSLINK wslp, array_meterp meterp, const int *          datap, int count));
WSDECL( int,   WSPutInteger64ArrayData,     ( WSLINK wslp, array_meterp meterp, const wsint64 *      datap, int count));
#elif WSINTERFACE == 3
WSDECL( int,   WSPutBinaryNumberArrayData,  ( WSLINK wslp, array_meterp meterp, void *         datap, long count, long type));
WSDECL( int,   WSPutByteArrayData,          ( WSLINK wslp, array_meterp meterp, unsigned char *datap, long count));

/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutShortIntegerArrayData  - WSPutInteger16ArrayData
WSPutIntegerArrayData       - WSPutInteger32ArrayData
WSPutLongIntegerArrayData   - WSPutInteger64ArrayData for 64-bit integers or WSPutInteger32ArrayData for 32-bit integers
*/
WSDECL( int,   WSPutShortIntegerArrayData,  ( WSLINK wslp, array_meterp meterp, short *        datap, long count));
WSDECL( int,   WSPutIntegerArrayData,       ( WSLINK wslp, array_meterp meterp, int *          datap, long count));
WSDECL( int,   WSPutLongIntegerArrayData,   ( WSLINK wslp, array_meterp meterp, long *         datap, long count));

WSDECL( int,   WSPutInteger16ArrayData,     ( WSLINK wslp, array_meterp meterp, short *        datap, int count));
WSDECL( int,   WSPutInteger32ArrayData,     ( WSLINK wslp, array_meterp meterp, int *          datap, int count));
WSDECL( int,   WSPutInteger64ArrayData,     ( WSLINK wslp, array_meterp meterp, wsint64 *      datap, int count));
#else
WSDECL( wsapi_result,   WSPutBinaryNumberArrayData,  ( WSLINK wslp, array_meterp meterp, voidp_ct     datap, long_st count, long type));
WSDECL( wsapi_result,   WSPutByteArrayData,          ( WSLINK wslp, array_meterp meterp, ucharp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSPutShortIntegerArrayData,  ( WSLINK wslp, array_meterp meterp, shortp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSPutIntegerArrayData,       ( WSLINK wslp, array_meterp meterp, intp_nt      datap, long_st count));
WSDECL( wsapi_result,   WSPutLongIntegerArrayData,   ( WSLINK wslp, array_meterp meterp, longp_nt     datap, long_st count));
#endif /* WSINTERFACE >= 3 */




#if WSINTERFACE >= 4
WSDECL( int,   WSPutFloatArrayData,         ( WSLINK wslp, array_meterp meterp, const float * datap, long count));
WSDECL( int,   WSPutDoubleArrayData,        ( WSLINK wslp, array_meterp meterp, const double *datap, long count));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSPutLongDoubleArrayData,    ( WSLINK wslp, array_meterp meterp, const wsextended_double *datap, long count));
#endif

WSDECL( int,   WSPutReal32ArrayData,        ( WSLINK wslp, array_meterp meterp, const float * datap, int count));
WSDECL( int,   WSPutReal64ArrayData,        ( WSLINK wslp, array_meterp meterp, const double *datap, int count));
WSDECL( int,   WSPutReal128ArrayData,       ( WSLINK wslp, array_meterp meterp, const wsextended_double *datap, int count));
#elif WSINTERFACE == 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutFloatArrayData      - WSPutReal32ArrayData
WSPutDoubleArrayData     - WSPutReal64ArrayData
WSPutLongDoubleArrayData - WSPutReal128ArrayData
*/
WSDECL( int,   WSPutFloatArrayData,         ( WSLINK wslp, array_meterp meterp, float * datap, long count));
WSDECL( int,   WSPutDoubleArrayData,        ( WSLINK wslp, array_meterp meterp, double *datap, long count));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSPutLongDoubleArrayData,    ( WSLINK wslp, array_meterp meterp, wsextended_double *datap, long count));
#endif

WSDECL( int,   WSPutReal32ArrayData,        ( WSLINK wslp, array_meterp meterp, float * datap, int count));
WSDECL( int,   WSPutReal64ArrayData,        ( WSLINK wslp, array_meterp meterp, double *datap, int count));
WSDECL( int,   WSPutReal128ArrayData,       ( WSLINK wslp, array_meterp meterp, wsextended_double *datap, int count));
#else
WSDECL( wsapi_result,   WSPutFloatArrayData,         ( WSLINK wslp, array_meterp meterp, floatp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSPutDoubleArrayData,        ( WSLINK wslp, array_meterp meterp, doublep_nt   datap, long_st count));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSPutLongDoubleArrayData,   ( WSLINK wslp, array_meterp meterp, extendedp_nt datap, long_st count));
#endif
#endif /* WSINTERFACE >= 4 */

#ifndef WS_USES_NEW_PUTBYTEARRAY_API
#define WS_USES_NEW_PUTBYTEARRAY_API 1
#endif


#if WSINTERFACE >= 4
WSDECL( int,   WSPutBinaryNumberArray,  ( WSLINK wslp, const void *         data, const long *dimp, const char **heads, long depth, long type));
WSDECL( int,   WSPutByteArray,          ( WSLINK wslp, const unsigned char *data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutShortIntegerArray,  ( WSLINK wslp, const short *        data, const long *dims, const char **heads, long depth));
WSDECL( int,   WSPutIntegerArray,       ( WSLINK wslp, const int *          data, const long *dims, const char **heads, long depth));
WSDECL( int,   WSPutLongIntegerArray,   ( WSLINK wslp, const long *         data, const long *dims, const char **heads, long depth));

WSDECL( int,   WSPutInteger8Array,      ( WSLINK wslp, const unsigned char *data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutInteger16Array,     ( WSLINK wslp, const short *        data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutInteger32Array,     ( WSLINK wslp, const int *          data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutInteger64Array,     ( WSLINK wslp, const wsint64 *      data, const int *dims, const char **heads, int depth));
#elif WSINTERFACE == 3
WSDECL( int,   WSPutBinaryNumberArray,  ( WSLINK wslp, void *         data, long *dimp, char **heads, long depth, long type));
WSDECL( int,   WSPutByteArray,          ( WSLINK wslp, unsigned char *data, int *dims, char **heads, int depth));

/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutShortIntegerArray  - WSPutInteger16Array
WSPutIntegerArray       - WSPutInteger32Array
WSPutLongIntegerArray   - WSPutInteger64Array for 64-bit integers or WSPutInteger32Array for 32-bit integers
*/
WSDECL( int,   WSPutShortIntegerArray,  ( WSLINK wslp, short *        data, long *dims, char **heads, long depth));
WSDECL( int,   WSPutIntegerArray,       ( WSLINK wslp, int *          data, long *dims, char **heads, long depth));
WSDECL( int,   WSPutLongIntegerArray,   ( WSLINK wslp, long *         data, long *dims, char **heads, long depth));

WSDECL( int,   WSPutInteger16Array,     ( WSLINK wslp, short *        data, int *dims, char **heads, int depth));
WSDECL( int,   WSPutInteger32Array,     ( WSLINK wslp, int *          data, int *dims, char **heads, int depth));
WSDECL( int,   WSPutInteger64Array,     ( WSLINK wslp, wsint64 *      data, int *dims, char **heads, int depth));
#else
WSDECL( wsapi_result,   WSPutBinaryNumberArray,  ( WSLINK wslp, voidp_ct     data, longp_st dimp, charpp_ct heads, long_st depth, long type));
WSDECL( wsapi_result,   WSPutByteArray,          ( WSLINK wslp, ucharp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( wsapi_result,   WSPutShortIntegerArray,  ( WSLINK wslp, shortp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( wsapi_result,   WSPutIntegerArray,       ( WSLINK wslp, intp_nt      data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( wsapi_result,   WSPutLongIntegerArray,   ( WSLINK wslp, longp_nt     data, longp_st dims, charpp_ct heads, long_st depth));
#endif /* WSINTERFACE >= 4 */


#if WSINTERFACE >= 4
WSDECL( int,   WSPutFloatArray,         ( WSLINK wslp, const float * data, const long *dims, const char **heads, long depth));
WSDECL( int,   WSPutDoubleArray,        ( WSLINK wslp, const double *data, const long *dims, const char **heads, long depth));
WSDECL( int,   WSPutRealArray,          ( WSLINK wslp, const double *data, const long *dims, const char **heads, long depth));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSPutLongDoubleArray,    ( WSLINK wslp, const wsextended_double *data, const long *dims, const char **heads, long depth));
#endif

WSDECL( int,   WSPutReal32Array,        ( WSLINK wslp, const float * data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutReal64Array,        ( WSLINK wslp, const double *data, const int *dims, const char **heads, int depth));
WSDECL( int,   WSPutReal128Array,       ( WSLINK wslp, const wsextended_double *data, const int *dims, const char **heads, int depth));
#elif WSINTERFACE == 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutFloatArray      - WSPutReal32Array
WSPutDoubleArray     - WSPutReal64Array
WSPutRealArray       - WSPutReal64Array
WSPutLongDoubleArray - WSPutReal128Array
*/
WSDECL( int,   WSPutFloatArray,         ( WSLINK wslp, float * data, long *dims, char **heads, long depth));
WSDECL( int,   WSPutDoubleArray,        ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
WSDECL( int,   WSPutRealArray,          ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSPutLongDoubleArray,    ( WSLINK wslp, wsextended_double *data, long *dims, char **heads, long depth));
#endif

WSDECL( int,   WSPutReal32Array,        ( WSLINK wslp, float * data, int *dims, char **heads, int depth));
WSDECL( int,   WSPutReal64Array,        ( WSLINK wslp, double *data, int *dims, char **heads, int depth));
WSDECL( int,   WSPutReal128Array,       ( WSLINK wslp, wsextended_double *data, int *dims, char **heads, int depth));
#else
WSDECL( wsapi_result,   WSPutFloatArray,         ( WSLINK wslp, floatp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( wsapi_result,   WSPutDoubleArray,        ( WSLINK wslp, doublep_nt   data, longp_st dims, charpp_ct heads, long_st depth));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSPutLongDoubleArray,   ( WSLINK wslp, extendedp_nt data, longp_st dims, charpp_ct heads, long_st depth));
#endif
#endif /* WSINTERFACE >= 4 */


#if WSINTERFACE >= 4
WSDECL( int,   WSPutBinaryNumberList, ( WSLINK wslp, const void *  data, long count, long type));
WSDECL( int,   WSPutIntegerList,      ( WSLINK wslp, const int *   data, long count));
WSDECL( int,   WSPutRealList,         ( WSLINK wslp, const double *data, long count));

WSDECL( int,   WSPutInteger8List,     ( WSLINK wslp, const unsigned char *data, int count));
WSDECL( int,   WSPutInteger16List,    ( WSLINK wslp, const short *   data, int count));
WSDECL( int,   WSPutInteger32List,    ( WSLINK wslp, const int *     data, int count));
WSDECL( int,   WSPutInteger64List,    ( WSLINK wslp, const wsint64 * data, int count));

WSDECL( int,   WSPutReal32List,       ( WSLINK wslp, const float * data, int count));
WSDECL( int,   WSPutReal64List,       ( WSLINK wslp, const double *data, int count));
WSDECL( int,   WSPutReal128List,      ( WSLINK wslp, const wsextended_double *data, int count));
#elif WSINTERFACE == 3
WSDECL( int,   WSPutBinaryNumberList, ( WSLINK wslp, void *  data, long count, long type));

/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSPutIntegerList  - WSPutInteger32List
WSPutRealList     - WSPutReal64List
*/
WSDECL( int,   WSPutIntegerList,      ( WSLINK wslp, int *   data, long count));
WSDECL( int,   WSPutRealList,         ( WSLINK wslp, double *data, long count));

WSDECL( int,   WSPutInteger16List,    ( WSLINK wslp, short *   data, int count));
WSDECL( int,   WSPutInteger32List,    ( WSLINK wslp, int *     data, int count));
WSDECL( int,   WSPutInteger64List,    ( WSLINK wslp, wsint64 * data, int count));

WSDECL( int,   WSPutReal32List,       ( WSLINK wslp, float * data, int count));
WSDECL( int,   WSPutReal64List,       ( WSLINK wslp, double *data, int count));
WSDECL( int,   WSPutReal128List,      ( WSLINK wslp, wsextended_double *data, int count));
#else
WSDECL( wsapi_result,   WSPutBinaryNumberList, ( WSLINK wslp, voidp_ct   data, long_st count, long type));
WSDECL( wsapi_result,   WSPutIntegerList,      ( WSLINK wslp, intp_nt    data, long_st count));
WSDECL( wsapi_result,   WSPutRealList,         ( WSLINK wslp, doublep_nt data, long_st count));
#endif /* WSINTERFACE >= 4 */


#if WSINTERFACE >= 3
#if WSINTERFACE >= 4
WSDECL( int, WSPutArrayType,             ( WSLINK wslp, WSLINK heads, long depth, array_meterpp meterpp));
WSDECL( int, WSReleasePutArrayState,     ( WSLINK wslp, WSLINK heads, array_meterp meterp));
#else
WSDECL( int, WSPutArrayType0,             ( WSLINK wslp, WSLINK heads, long depth, array_meterpp meterpp));
WSDECL( int, WSReleasePutArrayState0,     ( WSLINK wslp, WSLINK heads, array_meterp meterp));
#endif
#else
WSDECL( wsapi_result, WSPutArrayType0,             ( WSLINK wslp, WSLINK heads, long depth, array_meterpp meterpp));
WSDECL( wsapi_result, WSReleasePutArrayState0,     ( WSLINK wslp, WSLINK heads, array_meterp meterp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
#if WSINTERFACE >= 4
WSDECL( int, WSPutArrayLeaves,           ( WSLINK wslp, WSLINK heads, array_meterp meterp, WSLINK leaves, long count));
WSDECL( int, WSPutBinaryNumberArrayDataWithHeads, ( WSLINK wslp, WSLINK heads, array_meterp meterp, const void *datap, long count, long type));
#else
WSDECL( int, WSPutArrayLeaves0,           ( WSLINK wslp, WSLINK heads, array_meterp meterp, WSLINK leaves, long count));
WSDECL( int, WSPutBinaryNumberArrayData0, ( WSLINK wslp, WSLINK heads, array_meterp meterp, void *datap, long count, long type));
#endif
#else
#if WSINTERFACE > 1
WSDECL( wsapi_result, WSPutBinaryNumberArrayData0, ( WSLINK wslp, WSLINK heads, array_meterp meterp, voidp_ct datap, long_st count, long type));
WSDECL( wsapi_result, WSPutArrayLeaves0,           ( WSLINK wslp, WSLINK heads, array_meterp meterp, WSLINK leaves, long_st count));
#endif /* WSINTERFACE >= 1 */
#endif /* WSINTERFACE >= 3 */


WS_END_EXTERN_C




#ifndef _WSCAGET_H
#define _WSCAGET_H






#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter FAR * array_meterp;
typedef array_meterp FAR * array_meterpp;
#endif


#if WSINTERFACE < 3
#define WSGetRealArray    WSGetDoubleArray
#define WSDisownRealArray WSDisownDoubleArray
#endif

#endif /* _WSCAGET_H */



/* explicitly not protected by _WSCAGET_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( int,   WSGetArrayDimensions,       ( WSLINK wslp, array_meterp meterp));
WSDECL( int,   WSGetArrayType,             ( WSLINK wslp, array_meterp meterp));
#else
WSDECL( wsapi_result,   WSGetArrayDimensions,       ( WSLINK wslp, array_meterp meterp));
WSDECL( wsapi_token,    WSGetArrayType,             ( WSLINK wslp, array_meterp meterp));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSDECL( int,  WSGetBinaryNumberList, ( WSLINK wslp, void **datap, long *countp, long type));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetIntegerList  - WSGetInteger32List
WSGetRealList     - WSGetReal64List
*/
WSDECL( int,  WSGetIntegerList,      ( WSLINK wslp, int **datap, long *countp));
WSDECL( int,  WSGetRealList,         ( WSLINK wslp, double **datap, long *countp));

WSDECL( int,  WSGetInteger16List,      ( WSLINK wslp, short **   datap, int *countp));
WSDECL( int,  WSGetInteger32List,      ( WSLINK wslp, int **     datap, int *countp));
WSDECL( int,  WSGetInteger64List,      ( WSLINK wslp, wsint64 ** datap, int *countp));

WSDECL( int,  WSGetReal32List,         ( WSLINK wslp, float **                 datap, int *countp));
WSDECL( int,  WSGetReal64List,         ( WSLINK wslp, double **                datap, int *countp));
WSDECL( int,  WSGetReal128List,         ( WSLINK wslp, wsextended_double **    datap, int *countp));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSDisownBinaryNumberList  - WSReleaseBinaryNumberList
WSDisownIntegerList       - WSReleaseInteger32List
WSDisownRealList          - WSReleaseReal64List
*/

#if WSINTERFACE == 3
WSDECL( void, WSDisownBinaryNumberList, ( WSLINK wslp, void *data, long count, long type));
WSDECL( void, WSDisownIntegerList,      ( WSLINK wslp, int *data, long count));
WSDECL( void, WSDisownRealList,         ( WSLINK wslp, double *data, long count));
#endif

#if WSINTERFACE >= 4
WSDECL( void, WSReleaseIntegerList,     ( WSLINK wslp, int *data, long count));
WSDECL( void, WSReleaseRealList,        ( WSLINK wslp, double *data, long count));
#endif

WSDECL( void, WSReleaseBinaryNumberList,   ( WSLINK wslp, void *data, int count, long type));
WSDECL( void, WSReleaseInteger16List,      ( WSLINK wslp, short *data, int count));
WSDECL( void, WSReleaseInteger32List,      ( WSLINK wslp, int *data, int count));
WSDECL( void, WSReleaseInteger64List,      ( WSLINK wslp, wsint64 *data, int count));

WSDECL( void, WSReleaseReal32List,         ( WSLINK wslp, float *data, int count));
WSDECL( void, WSReleaseReal64List,         ( WSLINK wslp, double *data, int count));
WSDECL( void, WSReleaseReal128List,        ( WSLINK wslp, wsextended_double *data, int count));

WSDECL( int,   WSGetBinaryNumberArrayData,  ( WSLINK wslp, array_meterp meterp, void *datap, long count, long type));
WSDECL( int,   WSGetByteArrayData,          ( WSLINK wslp, array_meterp meterp, unsigned char * datap, long count));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetShortIntegerArrayData  - WSGetInteger16ArrayData
WSGetIntegerArrayData       - WSGetInteger32ArrayData
WSGetLongIntegerArrayData   - WSGetInteger64ArrayData for 64-bit integers or WSGetInteger32ArrayData for 32-bit integers.
*/

WSDECL( int,   WSGetShortIntegerArrayData,  ( WSLINK wslp, array_meterp meterp, short *         datap, long count));
WSDECL( int,   WSGetIntegerArrayData,       ( WSLINK wslp, array_meterp meterp, int *           datap, long count));
WSDECL( int,   WSGetLongIntegerArrayData,   ( WSLINK wslp, array_meterp meterp, long *          datap, long count));

WSDECL( int,   WSGetInteger16ArrayData,     ( WSLINK wslp, array_meterp meterp, short *         datap, int count));
WSDECL( int,   WSGetInteger32ArrayData,     ( WSLINK wslp, array_meterp meterp, int *           datap, int count));
WSDECL( int,   WSGetInteger64ArrayData,     ( WSLINK wslp, array_meterp meterp, wsint64 *       datap, int count));
#else
WSDECL( wsapi_result,  WSGetBinaryNumberList, ( WSLINK wslp, voidpp_ct   datap, longp_st countp, long type));
WSDECL( wsapi_result,  WSGetIntegerList,      ( WSLINK wslp, intpp_nt    datap, longp_st countp));
WSDECL( wsapi_result,  WSGetRealList,         ( WSLINK wslp, doublepp_nt datap, longp_st countp));

WSDECL( void, WSDisownBinaryNumberList, ( WSLINK wslp, voidp_ct   data, long_st count, long type));
WSDECL( void, WSDisownIntegerList,      ( WSLINK wslp, intp_nt    data, long_st count));
WSDECL( void, WSDisownRealList,         ( WSLINK wslp, doublep_nt data, long_st count));

WSDECL( wsapi_result,   WSGetBinaryNumberArrayData,  ( WSLINK wslp, array_meterp meterp, voidp_ct     datap, long_st count, long type));
WSDECL( wsapi_result,   WSGetByteArrayData,          ( WSLINK wslp, array_meterp meterp, ucharp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSGetShortIntegerArrayData,  ( WSLINK wslp, array_meterp meterp, shortp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSGetIntegerArrayData,       ( WSLINK wslp, array_meterp meterp, intp_nt      datap, long_st count));
WSDECL( wsapi_result,   WSGetLongIntegerArrayData,   ( WSLINK wslp, array_meterp meterp, longp_nt     datap, long_st count));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetFloatArrayData        - WSGetReal32ArrayData
WSGetDoubleArrayData       - WSGetReal64ArrayData
WSGetLongDoubleArrayData   - WSGetReal128ArrayData
*/

WSDECL( int,   WSGetFloatArrayData,         ( WSLINK wslp, array_meterp meterp, float *datap, long count));
WSDECL( int,   WSGetDoubleArrayData,        ( WSLINK wslp, array_meterp meterp, double *datap, long count));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSGetLongDoubleArrayData,   ( WSLINK wslp, array_meterp meterp, wsextended_double *datap, long count));
#endif

WSDECL( int,   WSGetReal32ArrayData,         ( WSLINK wslp, array_meterp meterp, float *datap, int count));
WSDECL( int,   WSGetReal64ArrayData,        ( WSLINK wslp, array_meterp meterp, double *datap, int count));
WSDECL( int,   WSGetReal128ArrayData,   ( WSLINK wslp, array_meterp meterp, wsextended_double *datap, int count));
#else
WSDECL( wsapi_result,   WSGetFloatArrayData,         ( WSLINK wslp, array_meterp meterp, floatp_nt    datap, long_st count));
WSDECL( wsapi_result,   WSGetDoubleArrayData,        ( WSLINK wslp, array_meterp meterp, doublep_nt   datap, long_st count));

#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSGetLongDoubleArrayData,   ( WSLINK wslp, array_meterp meterp, extendedp_nt datap, long_st count));
#endif

#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL(int, WSGetInteger8List, (WSLINK wslp, unsigned char **datap, int *countp));
WSDECL(int, WSGetInteger8ArrayData, (WSLINK wslp, array_meterp meterp, unsigned char *datap, int count));
WSDECL(void, WSReleaseInteger8List, (WSLINK wslp, unsigned char *data, int count));
#endif

#if WSINTERFACE >= 3

#if WSINTERFACE == 3
WSDECL( int,   WSGetArrayType0,             ( WSLINK wslp, WSLINK heads, array_meterpp meterpp, long *depthp, wsapi__token *leaf_tokp));
#else /* WSINTERFACE >= 4 */
WSDECL( int,   WSGetArrayTypeWithDepthAndLeafType, ( WSLINK wslp, WSLINK heads, array_meterpp meterpp,
	long *depthp, wsapi__token *leaf_tokp));
#endif

#if WSINTERFACE >= 4
WSDECL( int,   WSGetBinaryNumberArrayDataWithHeads, ( WSLINK wslp, WSLINK heads, array_meterp  meterp, void *datap, long *countp, long type));
#else
WSDECL( int,   WSGetBinaryNumberArrayData0, ( WSLINK wslp, WSLINK heads, array_meterp  meterp, void *datap, long *countp, long type));
#endif

#if WSINTERFACE >= 4
WSDECL( void,  WSReleaseGetArrayState,     ( WSLINK wslp, WSLINK heads, array_meterp  meterp));
#else
WSDECL( void,  WSReleaseGetArrayState0,     ( WSLINK wslp, WSLINK heads, array_meterp  meterp));
#endif

#if WSINTERFACE == 3
WSDECL( int,   WSGetBinaryNumberArray0,   ( WSLINK wslp, void **datap, long **dimpp, char ***headsp, long *depthp, long type, wsapi__token *leaf_tokp));
#else /* WSINTERFACE >= 4 */
WSDECL( int,   WSGetBinaryNumberArrayWithLeafType,   ( WSLINK wslp, void **datap, long **dimpp, char ***headsp, long *depthp, long type, wsapi__token *leaf_tokp));
#endif

#else
#if WSINTERFACE > 1
WSDECL( wsapi_result,   WSGetArrayType0,             ( WSLINK wslp, WSLINK heads, array_meterpp meterpp, longp_st depthp, wsapi__tokenp leaf_tokp));
WSDECL( wsapi_result,   WSGetBinaryNumberArrayData0, ( WSLINK wslp, WSLINK heads, array_meterp  meterp, voidp_ct datap, longp_st countp, long type));
WSDECL( void,           WSReleaseGetArrayState0,     ( WSLINK wslp, WSLINK heads, array_meterp  meterp));

WSDECL( wsapi_result,   WSGetBinaryNumberArray0,   ( WSLINK wslp, voidpp_ct     datap, longpp_st dimpp, charppp_ct headsp, longp_st depthp, long type, wsapi__tokenp leaf_tokp));
#endif /* WSINTERFACE > 1 */
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE >= 3
WSDECL( int,   WSGetBinaryNumberArray,    ( WSLINK wslp, void **          datap, long **dimpp, char ***headsp, long *depthp, long type));
WSDECL( int,   WSGetByteArray,            ( WSLINK wslp, unsigned char ** datap, int **dimsp, char ***headsp, int *depthp));

/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetShortIntegerArray   - WSGetInteger16Array
WSGetIntegerArray        - WSGetInteger32Array
WSGetLongIntegerArray    - WSGetInteger64Array for 64-bit integers or WSGetInteger32Array from 32-bit integers.
*/

WSDECL( int,   WSGetShortIntegerArray,    ( WSLINK wslp, short **         datap, long **dimsp, char ***headsp, long *depthp));
WSDECL( int,   WSGetIntegerArray,         ( WSLINK wslp, int **           datap, long **dimsp, char ***headsp, long *depthp));
WSDECL( int,   WSGetLongIntegerArray,     ( WSLINK wslp, long **          datap, long **dimsp, char ***headsp, long *depthp));

WSDECL( int,   WSGetInteger16Array,       ( WSLINK wslp, short **         datap, int **dimsp, char ***headsp, int *depthp));
WSDECL( int,   WSGetInteger32Array,       ( WSLINK wslp, int **           datap, int **dimsp, char ***headsp, int *depthp));
WSDECL( int,   WSGetInteger64Array,       ( WSLINK wslp, wsint64 **       datap, int **dimsp, char ***headsp, int *depthp));
#else
WSDECL( wsapi_result,   WSGetBinaryNumberArray,    ( WSLINK wslp, voidpp_ct     datap, longpp_st dimpp, charppp_ct headsp, longp_st depthp, long type));
WSDECL( wsapi_result,   WSGetByteArray,            ( WSLINK wslp, ucharpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
WSDECL( wsapi_result,   WSGetShortIntegerArray,    ( WSLINK wslp, shortpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
WSDECL( wsapi_result,   WSGetIntegerArray,         ( WSLINK wslp, intpp_nt      datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
WSDECL( wsapi_result,   WSGetLongIntegerArray,     ( WSLINK wslp, longpp_nt     datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE >= 4
WSDECL(int,  WSGetInteger8Array,  (WSLINK wslp, unsigned char **datap, int **dimsp, char ***headsp, int *depthp));
#endif

#if WSINTERFACE >= 3
/*
As of WSINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

WSGetFloatArray      - WSGetReal32Array
WSGetDoubleArray     - WSGetReal64Array
WSGetRealArray       - WSGetReal64Array
WSGetLongDoubleArray - WSGetReal128Array
*/

WSDECL( int,   WSGetFloatArray,           ( WSLINK wslp, float ** datap, long **dimsp, char ***headsp, long *depthp));
WSDECL( int,   WSGetDoubleArray,          ( WSLINK wslp, double **datap, long **dimsp, char ***headsp, long *depthp));
WSDECL( int,   WSGetRealArray,            ( WSLINK wslp, double **datap, long **dimsp, char ***headsp, long *depthp));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( int,   WSGetLongDoubleArray,      ( WSLINK wslp, wsextended_double **datap, long **dimsp, char ***headsp, long *depthp));
#endif

WSDECL( int,   WSGetReal32Array,          ( WSLINK wslp, float ** datap, int **dimsp, char ***headsp, int *depthp));
WSDECL( int,   WSGetReal64Array,          ( WSLINK wslp, double **datap, int **dimsp, char ***headsp, int *depthp));
WSDECL( int,   WSGetReal128Array,         ( WSLINK wslp, wsextended_double **datap, int **dimsp, char ***headsp, int *depthp));
#else
WSDECL( wsapi_result,   WSGetDoubleArray,          ( WSLINK wslp, doublepp_nt   datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
WSDECL( wsapi_result,   WSGetFloatArray,           ( WSLINK wslp, floatpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
#if CC_SUPPORTS_LONG_DOUBLE
WSDECL( wsapi_result,   WSGetLongDoubleArray,      ( WSLINK wslp, extendedpp_nt datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
#endif
#endif /* WSINTERFACE >= 3*/


#if WSINTERFACE >= 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSDisownBinaryNumberArray - WSReleaseBinaryNumberArray
WSDisownByteArray         - WSReleaseByteArray
WSDisownShortIntegerArray - WSReleaseInteger16Array
WSDisownIntegerArray      - WSReleaseInteger32Array
WSDisownLongIntegerArray  - Use WSReleaseInteger32Array for 32-bit integers or WSReleaseInteger64Array for 64-bit integers.
*/

#if WSINTERFACE == 3
WSDECL( void,           WSDisownBinaryNumberArray,  ( WSLINK wslp, void *         data, long *dimp, char **heads, long len, long type));
#endif

#if WSINTERFACE == 3
WSDECL( void,           WSDisownByteArray,          ( WSLINK wslp, unsigned char *data, long *dims, char **heads, long depth));
WSDECL( void,           WSDisownShortIntegerArray,  ( WSLINK wslp, short *        data, long *dims, char **heads, long depth));
WSDECL( void,           WSDisownIntegerArray,       ( WSLINK wslp, int *          data, long *dims, char **heads, long depth));
WSDECL( void,           WSDisownLongIntegerArray,   ( WSLINK wslp, long *         data, long *dims, char **heads, long depth));
#endif

#if WSINTERFACE >= 4
WSDECL( void,           WSReleaseShortIntegerArray, ( WSLINK wslp, short *        data, long *dims, char **heads, long depth));
WSDECL( void,           WSReleaseIntegerArray,      ( WSLINK wslp, int *          data, long *dims, char **heads, long depth));
WSDECL( void,           WSReleaseLongIntegerArray,  ( WSLINK wslp, long *         data, long *dims, char **heads, long depth));
#endif

WSDECL( void,           WSReleaseBinaryNumberArray,  ( WSLINK wslp, void *         data, int *dimp, char **heads, int len, long type));
WSDECL( void,           WSReleaseByteArray,          ( WSLINK wslp, unsigned char *data, int *dims, char **heads, int depth));
WSDECL( void,           WSReleaseInteger16Array,     ( WSLINK wslp, short *        data, int *dims, char **heads, int depth));
WSDECL( void,           WSReleaseInteger32Array,     ( WSLINK wslp, int *          data, int *dims, char **heads, int depth));
WSDECL( void,           WSReleaseInteger64Array,     ( WSLINK wslp, wsint64 *      data, int *dims, char **heads, int depth));
#else
WSDECL( void,           WSDisownBinaryNumberArray,  ( WSLINK wslp, voidp_ct     data, longp_st dimp, charpp_ct heads, long_st len, long type));
WSDECL( void,           WSDisownByteArray,          ( WSLINK wslp, ucharp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( void,           WSDisownShortIntegerArray,  ( WSLINK wslp, shortp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( void,           WSDisownIntegerArray,       ( WSLINK wslp, intp_nt      data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( void,           WSDisownLongIntegerArray,   ( WSLINK wslp, longp_nt     data, longp_st dims, charpp_ct heads, long_st depth));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WSDECL(void,   WSReleaseInteger8Array,   (WSLINK wslp, unsigned char *data, int *dimp, char **heads, int depth));
#endif

#if WSINTERFACE >= 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSDisownFloatArray  - WSReleaseReal32Array
WSDisownDoubleArray - WSReleaseReal64Array
WSDisownRealArray   - WSReleaseReal64Array
*/

#if WSINTERFACE == 3
WSDECL( void,           WSDisownFloatArray,         ( WSLINK wslp, float * data, long *dims, char **heads, long depth));
WSDECL( void,           WSDisownDoubleArray,        ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
WSDECL( void,           WSDisownRealArray,          ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
#endif

#if WSINTERFACE >= 4
WSDECL( void,           WSReleaseFloatArray,         ( WSLINK wslp, float * data, long *dims, char **heads, long depth));
WSDECL( void,           WSReleaseDoubleArray,        ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
WSDECL( void,           WSReleaseRealArray,          ( WSLINK wslp, double *data, long *dims, char **heads, long depth));
#endif

WSDECL( void,           WSReleaseReal32Array,         ( WSLINK wslp, float * data, int *dims, char **heads, int depth));
WSDECL( void,           WSReleaseReal64Array,          ( WSLINK wslp, double *data, int *dims, char **heads, int depth));
#else
WSDECL( void,           WSDisownFloatArray,         ( WSLINK wslp, floatp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
WSDECL( void,           WSDisownDoubleArray,        ( WSLINK wslp, doublep_nt   data, longp_st dims, charpp_ct heads, long_st depth));
#endif /* WSINTERFACE >= 3 */


#if CC_SUPPORTS_LONG_DOUBLE
#if WSINTERFACE >= 3
/* 
As of WSINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

WSDisownLongDoubleArray - WSReleaseReal128Array
*/

#if WSINTERFACE == 3
WSDECL( void,           WSDisownLongDoubleArray,   ( WSLINK wslp, wsextended_double *data, long *dims, char **heads, long depth));
#endif

WSDECL( void,           WSReleaseReal128Array,     ( WSLINK wslp, wsextended_double *data, int *dims, char **heads, int depth));

#if WSINTERFACE >= 4
WSDECL( void,           WSReleaseLongDoubleArray,  ( WSLINK wslp, wsextended_double *data, long *dims, char **heads, long depth));
#endif

#else
WSDECL( void,           WSDisownLongDoubleArray,   ( WSLINK wslp, extendedp_nt data, longp_st dims, charpp_ct heads, long_st depth));
#endif /* WSINTERFACE >= 3 */
#endif /* CC_SUPPORTS_LONG_DOUBLE */


WS_END_EXTERN_C



/*************** Unicode Container for mprep templates ***************/


#ifndef WSUNICODECONTAINER_HPP
#define WSUNICODECONTAINER_HPP







WS_EXTERN_C

#if WSINTERFACE >= 4
enum WSUnicodeContainerType
{
	UCS2ContainerType,
	UTF8ContainerType,
	UTF16ContainerType,
	UTF32ContainerType
};

typedef struct _WSUnicodeContainer
{
	union _pointer
	{
		unsigned short *ucs2;
		unsigned char *utf8;
		unsigned short *utf16;
		unsigned int *utf32;
	} pointer;

	int length;
	enum WSUnicodeContainerType type;
} WSUnicodeContainer;

#define WSUCS2String(container)   container->pointer.ucs2
#define WSUTF8String(container)   container->pointer.utf8
#define WSUTF16String(container)  container->pointer.utf16
#define WSUTF32String(container)  container->pointer.utf32
#define WSUnicodeStringLength(container) container->length
#define WSUnicodeStringType(container) container->type

WSDECL(WSUnicodeContainer *, WSNewUnicodeContainer, (const void *string, int length, enum WSUnicodeContainerType type));
WSDECL(void,    WSReleaseUnicodeContainer, (WSUnicodeContainer *string));

#endif /* WSINTERFACE >= 4 */

WS_END_EXTERN_C


#endif /* WSUNICODECONTAINER_HPP */




/*************** seeking, transfering  and synchronization ***************/

#ifndef _WSMARK_H
#define _WSMARK_H




#endif /* _WSMARK_H */

/* explicitly not protected by _WSMARK_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C
WSDECL( WSINKMark,  WSCreateMark,  ( WSLINK wslp));
#if WSINTERFACE >= 3
WSDECL( WSINKMark,  WSSeekToMark,  ( WSLINK wslp, WSINKMark mark, int index));
WSDECL( WSINKMark,  WSSeekMark,    ( WSLINK wslp, WSINKMark mark, int index));
#else
WSDECL( WSINKMark,  WSSeekToMark,  ( WSLINK wslp, WSINKMark mark, long index));
WSDECL( WSINKMark,  WSSeekMark,    ( WSLINK wslp, WSINKMark mark, long index));
#endif /* WSINTERFACE >= 3 */
WSDECL( void,       WSDestroyMark, ( WSLINK wslp, WSINKMark mark));
WS_END_EXTERN_C




#ifndef _WSXFER_H
#define _WSXFER_H





#endif /* _WSXFER_H */

/* explicitly not protected by _WSXFER_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C

#ifndef WSINTERFACE
/*\ syntax\ error\ */\ )
#endif

#if WSINTERFACE >= 3
WSDECL( int, WSTransferExpression, ( WSLINK dwslp, WSLINK swslp));
WSDECL( int, WSTransferToEndOfLoopbackLink, ( WSLINK dwslp, WSLINK swslp));
#else
WSDECL( wsapi_result, WSTransferExpression, ( WSLINK dwslp, WSLINK swslp));
WSDECL( wsapi_result, WSTransferToEndOfLoopbackLink, ( WSLINK dwslp, WSLINK swslp));
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE > 1
#if WSINTERFACE <= 3
#if WSINTERFACE == 3
WSDECL( int, WSTransfer0, ( WSLINK dwslp, WSLINK swslp, unsigned long sequence_no));
#else
WSDECL( wsapi_result, WSTransfer0, ( WSLINK dwslp, WSLINK swslp, ulong_ct sequence_no));
#endif
#endif /* WSINTERFACE <= 3 */
#else
static wsapi_result WSTransfer0 (WSLINK dwslp, WSLINK swslp, ulong_ct sequence_no);
#endif /* WSINTERFACE > 1 */
WS_END_EXTERN_C




#ifndef _WSSYNC_H
#define _WSSYNC_H





#endif /* _WSSYNC_H */

/* explicitly not protected by _WSSYNC_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C
/* in response to a reset message */
#if WSINTERFACE >= 3
WSDECL( int, WSForwardReset, ( WSLINK wslp, unsigned long marker));
#else
WSDECL( wsapi_result, WSForwardReset, ( WSLINK wslp, ulong_ct marker));
#endif /* WSINTERFACE >= 3 */


#if WSINTERFACE >= 3
WSDECL( int, WSAlign,        ( WSLINK lwslp, WSLINK rwslp));
#else
WSDECL( wsapi_result, WSAlign,        ( WSLINK lwslp, WSLINK rwslp));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C



/*************************************************************/


#ifndef _WSPKT_H
#define _WSPKT_H

/*************** Mathematica packet interface ***************/

			/* WSNextPacket returns one of... */



#ifndef _WSPKTNO_H
#define _WSPKTNO_H

#define ILLEGALPKT      0

#define CALLPKT         7
#define EVALUATEPKT    13
#define RETURNPKT       3

#define INPUTNAMEPKT    8
#define ENTERTEXTPKT   14
#define ENTEREXPRPKT   15
#define OUTPUTNAMEPKT   9
#define RETURNTEXTPKT   4
#define RETURNEXPRPKT  16

#define DISPLAYPKT     11
#define DISPLAYENDPKT  12

#define MESSAGEPKT      5
#define TEXTPKT         2

#define INPUTPKT        1
#define INPUTSTRPKT    21
#define MENUPKT         6
#define SYNTAXPKT      10

#define SUSPENDPKT     17
#define RESUMEPKT      18

#define BEGINDLGPKT    19
#define ENDDLGPKT      20

#define FIRSTUSERPKT  128
#define LASTUSERPKT   255


#endif /* _WSPKTNO_H */




#endif /* _WSPKT_H */

/* explicitly not protected by _WSPKT_H in case WSDECL is redefined for multiple inclusion */

WS_EXTERN_C
#if WSINTERFACE >= 3
WSDECL( int,  WSNextPacket, ( WSLINK wslp));
#else
WSDECL( wsapi_packet,  WSNextPacket, ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C




#ifndef _WSALERT_H
#define _WSALERT_H







WS_EXTERN_C
/*************** User interaction--for internal use only ***************/

#if WIN64_WSTP
typedef __int64 wsdlg_result;
#else
typedef long wsdlg_result;
#endif

#if WSINTERFACE >= 3
WSDPROC( wsdlg_result, WSAlertProcPtr,             ( WSEnvironment env, const char *message));
WSDPROC( wsdlg_result, WSRequestProcPtr,           ( WSEnvironment env, const char *prompt, char *response, long sizeof_response));
WSDPROC( wsdlg_result, WSConfirmProcPtr,           ( WSEnvironment env, const char *question, wsdlg_result default_answer));
WSDPROC( wsdlg_result, WSRequestArgvProcPtr,       ( WSEnvironment env, char **argv, long cardof_argv, char *buf, long sizeof_buf));
#else
WSDPROC( wsdlg_result, WSAlertProcPtr,             ( WSEnvironment env, kcharp_ct message));
WSDPROC( wsdlg_result, WSRequestProcPtr,           ( WSEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
WSDPROC( wsdlg_result, WSConfirmProcPtr,           ( WSEnvironment env, kcharp_ct question, wsdlg_result default_answer));
WSDPROC( wsdlg_result, WSRequestArgvProcPtr,       ( WSEnvironment env, charpp_ct argv, long cardof_argv, charp_ct buf, long sizeof_buf));
#endif /* WSINTERFACE >= 3 */
WSDPROC( wsdlg_result, WSRequestToInteractProcPtr, ( WSEnvironment env, wsdlg_result wait_for_permission));
WSDPROC( wsdlg_result, WSDialogProcPtr,            ( WSEnvironment env));


typedef WSDialogProcPtr WSDialogUPP;
typedef WSAlertProcPtr WSAlertUPP;
typedef WSRequestProcPtr WSRequestUPP;
typedef WSConfirmProcPtr WSConfirmUPP;
typedef WSRequestArgvProcPtr WSRequestArgvUPP;
typedef WSRequestToInteractProcPtr WSRequestToInteractUPP;
#define NewWSAlertProc(userRoutine) WSAlertCast((userRoutine))
#define NewWSRequestProc(userRoutine) WSRequestCast((userRoutine))
#define NewWSConfirmProc(userRoutine) WSConfirmCast((userRoutine))
#define NewWSRequestArgvProc(userRoutine) WSRequestArgvCast((userRoutine))
#define NewWSRequestToInteractProc(userRoutine) WSRequestToInteractCast((userRoutine))

typedef WSAlertUPP WSAlertFunctionType;
typedef WSRequestUPP WSRequestFunctionType;
typedef WSConfirmUPP WSConfirmFunctionType;
typedef WSRequestArgvUPP WSRequestArgvFunctionType;
typedef WSRequestToInteractUPP WSRequestToInteractFunctionType;
typedef WSDialogUPP WSDialogFunctionType;



/* 
	WSDDECL( wsdlg_result, alert_user, ( WSEnvironment env, kcharp_ct message));
	WSDDEFN( wsdlg_result, alert_user, ( WSEnvironment env, kcharp_ct message))
	{
		fprintf( stderr, "%s\n", message);
	}


	...
	WSDialogFunctionType f = NewWSAlertProc(alert_user);
	WSSetDialogFunction( ep, WSAlertFunction, f);
	...
	or
	...
	WSSetDialogFunction( ep, WSAlertFunction, NewWSAlertProc(alert_user));
	...
*/



enum {	WSAlertFunction = 1, WSRequestFunction, WSConfirmFunction,
	WSRequestArgvFunction, WSRequestToInteractFunction };


#define WS_DEFAULT_DIALOG ( (WSDialogFunctionType) 1)
#define WS_IGNORE_DIALOG ( (WSDialogFunctionType) 0)
#define WS_SUPPRESS_DIALOG WS_IGNORE_DIALOG



#if WINDOWS_WSTP

#ifndef _WSWIN_H
#define _WSWIN_H





WS_EXTERN_C
#if WSINTERFACE >= 3
WSDDECL( wsdlg_result, WSAlert_win,   ( WSEnvironment ep, const char *alertstr));
WSDDECL( wsdlg_result, WSRequest_win, ( WSEnvironment ep, const char *prompt, char *response, long n));
WSDDECL( wsdlg_result, WSConfirm_win, ( WSEnvironment ep, const char *okcancelquest, wsdlg_result default_answer));
#else
WSDDECL( wsdlg_result, WSAlert_win,   ( WSEnvironment ep, kcharp_ct alertstr));
WSDDECL( wsdlg_result, WSRequest_win, ( WSEnvironment ep, kcharp_ct prompt, charp_ct response, long n));
WSDDECL( wsdlg_result, WSConfirm_win, ( WSEnvironment ep, kcharp_ct okcancelquest, wsdlg_result default_answer));
#endif /* WSINTERFACE >= 3 */
WSDDECL( wsdlg_result, WSPermit_win,  ( WSEnvironment ep, wsdlg_result wait));
WS_END_EXTERN_C


#define DLG_LINKNAME                101
#define DLG_TEXT                    102
#define RIDOK                       1
#define RIDCANCEL                   104

#endif /* _WSWIN_H */


#define WSALERT         WSAlert_win
#define WSREQUEST       WSRequest_win
#define WSCONFIRM       WSConfirm_win
#define WSPERMIT        WSPermit_win
#define WSREQUESTARGV	default_request_argv
#endif

#if UNIX_WSTP
#if DARWIN_WSTP  && ! defined(IOS_WSTP) && defined (USE_CF_DIALOGS)

#ifndef _WSDARWIN_H
#define _WSDARWIN_H





WS_EXTERN_C

#if WSINTERFACE >= 3
WSDDECL( wsdlg_result, WSAlert_darwin,   ( WSEnvironment env, const char *message));
WSDDECL( wsdlg_result, WSRequest_darwin, ( WSEnvironment env, const char *prompt, char *response, long sizeof_response));
WSDDECL( wsdlg_result, WSConfirm_darwin, ( WSEnvironment env, const char *question, wsdlg_result default_answer));
#else
WSDDECL( wsdlg_result, WSAlert_darwin,   ( WSEnvironment env, kcharp_ct message));
WSDDECL( wsdlg_result, WSRequest_darwin, ( WSEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
WSDDECL( wsdlg_result, WSConfirm_darwin, ( WSEnvironment env, kcharp_ct question, wsdlg_result default_answer));
#endif /* WSINTERFACE >= 3 */
WSDDECL( wsdlg_result, WSPermit_darwin,  ( WSEnvironment env, wsdlg_result wait_for_permission));
WSDDECL( wsdlg_result, WSDontPermit_darwin, ( WSEnvironment ep, wsdlg_result wait_for_permission));

WS_END_EXTERN_C

#endif /* _WSDARWIN_H */


#define WSALERT  	WSAlert_darwin
#define WSREQUEST	WSRequest_darwin
#define WSCONFIRM	WSConfirm_darwin
#define WSPERMIT 	WSPermit_darwin
#define WSREQUESTARGV	default_request_argv
#else /* !(DARWIN_WSTP && ! defined(IOS_WSTP) && defined (USE_CF_DIALOGS)) */

#ifndef _WSUNIX_H
#define _WSUNIX_H





WS_EXTERN_C

#if WSINTERFACE >= 3
WSDDECL( wsdlg_result, WSAlert_unix,   ( WSEnvironment env, const char *message));
WSDDECL( wsdlg_result, WSRequest_unix, ( WSEnvironment env, const char *prompt, char *response, long sizeof_response));
WSDDECL( wsdlg_result, WSConfirm_unix, ( WSEnvironment env, const char *question, wsdlg_result default_answer));
#else
WSDDECL( wsdlg_result, WSAlert_unix,   ( WSEnvironment env, kcharp_ct message));
WSDDECL( wsdlg_result, WSRequest_unix, ( WSEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
WSDDECL( wsdlg_result, WSConfirm_unix, ( WSEnvironment env, kcharp_ct question, wsdlg_result default_answer));
#endif /* WSINTERFACE >= 3 */
WSDDECL( wsdlg_result, WSPermit_unix,  ( WSEnvironment env, wsdlg_result wait_for_permission));

WS_END_EXTERN_C

#endif /* _WSUNIX_H */


#define WSALERT  	WSAlert_unix
#define WSREQUEST	WSRequest_unix
#define WSCONFIRM	WSConfirm_unix
#define WSPERMIT 	WSPermit_unix
#define WSREQUESTARGV	default_request_argv
#endif /* DARWIN_WSTP && ! defined(IOS_WSTP) && defined (USE_CF_DIALOGS) */
#endif


#if WSINTERFACE >= 3
WSDDECL( wsdlg_result, default_request_argv, ( WSEnvironment ep, char **argv, long len, char *buff, long size));
#else
WSDDECL( wsdlg_result, default_request_argv, ( WSEnvironment ep, charpp_ct argv, long len, charp_ct buff, long size));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C

#endif /* _WSALERT_H */

WS_EXTERN_C

#if WSINTERFACE >= 3
WSDECL( wsdlg_result,  WSAlert,             ( WSEnvironment env, const char *message));
WSDECL( wsdlg_result,  WSRequest,           ( WSEnvironment env, const char *prompt, char *response, long sizeof_response)); /* initialize response with default*/
WSDECL( wsdlg_result,  WSConfirm,           ( WSEnvironment env, const char *question, wsdlg_result default_answer));
WSDECL( wsdlg_result,  WSRequestArgv,       ( WSEnvironment env, char **argv, long cardof_argv, char *buff, long size));
#else
WSDECL( wsdlg_result,  WSAlert,             ( WSEnvironment env, kcharp_ct message));
WSDECL( wsdlg_result,  WSRequest,           ( WSEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response)); /* initialize response with default*/
WSDECL( wsdlg_result,  WSConfirm,           ( WSEnvironment env, kcharp_ct question, wsdlg_result default_answer));
WSDECL( wsdlg_result,  WSRequestArgv,       ( WSEnvironment env, charpp_ct argv, long cardof_argv, charp_ct buff, long size));
#endif /* WSINTERFACE >= 3 */

WSDECL( wsdlg_result,  WSRequestToInteract, ( WSEnvironment env, wsdlg_result wait_for_permission));
#if WSINTERFACE >= 3
WSDECL( int,  WSSetDialogFunction, ( WSEnvironment env, long funcnum, WSDialogFunctionType func));
#else
WSDECL( wsapi_result,  WSSetDialogFunction, ( WSEnvironment env, long funcnum, WSDialogFunctionType func));
#endif /* WSINTERFACE >= 3 */

/* just some type-safe casts */
WSDECL( WSDialogProcPtr, WSAlertCast, ( WSAlertProcPtr f));
WSDECL( WSDialogProcPtr, WSRequestCast, ( WSRequestProcPtr f));
WSDECL( WSDialogProcPtr, WSConfirmCast, ( WSConfirmProcPtr f));
WSDECL( WSDialogProcPtr, WSRequestArgvCast, ( WSRequestArgvProcPtr f));
WSDECL( WSDialogProcPtr, WSRequestToInteractCast, ( WSRequestToInteractProcPtr f));
WS_END_EXTERN_C



/*************************************************************/


#ifndef _WSREADY_H
#define _WSREADY_H

#ifndef _WSTIME_H
#define _WSTIME_H

typedef struct _wstimeval{
	unsigned long tv_sec;
	unsigned long tv_usec;
} wstimeval;


#endif /* _WSTIME_H */




WS_EXTERN_C
#if WSINTERFACE >= 3
WSDECL( int,   WSReady,            ( WSLINK wslp));
#else
WSDECL( wsapi_result,   WSReady,            ( WSLINK wslp));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C

#if WSINTERFACE >= 3


#if defined(__cplusplus)
#define WSINFINITEWAIT static_cast<wstimeval>({-1, -1})
#else
#define WSINFINITEWAIT {-1, -1}
#endif /* __cplusplus */

#define WSREADYPARALLELERROR -1
#define WSREADYPARALLELTIMEDOUT -2
#define WSREADYPARALLELINVALIDARGUMENT -3

typedef void *WSREADYPARALLELENV;

WS_EXTERN_C
WSDECL(int,              WSReadyParallel, (WSENV, WSLINK *, int, wstimeval));
WS_END_EXTERN_C

#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 4
WS_EXTERN_C
WSWPROC(int, WSLinkWaitCallBackObject, (WSLINK, void *));

WSDECL(int, WSWaitForLinkActivity, (WSLINK wslp));
WSDECL(int, WSWaitForLinkActivityWithCallback, (WSLINK wslp, WSLinkWaitCallBackObject callback));
WS_END_EXTERN_C


#define WSWAITSUCCESS         1
#define WSWAITERROR           2
#define WSWAITTIMEOUT         3
#define WSWAITCALLBACKABORTED 4
#endif

#endif /* _WSREADY_H */



/********************************************************/


#ifndef _WSTM_H
#define _WSTM_H




/*************** Template interface ***************/

/* The following are useful only when using template files as
 * their definitions are produced by mprep.
 */

WS_EXTERN_C
extern WSLINK stdlink;
extern WSEnvironment stdenv;

extern WSYieldFunctionObject stdyielder;
extern WSMessageHandlerObject stdhandler;

#if WSINTERFACE >= 3
extern int WSMain P((int, char **)); /* pass in argc and argv */
extern int WSMainString P(( char *commandline));
#else
extern int WSMain P((int, charpp_ct)); /* pass in argc and argv */
extern int WSMainString P(( charp_ct commandline));
#endif /* WSINTERFACE >= 3 */

extern int WSMainArgv P(( char** argv, char** argv_end)); /* note not FAR pointers */
            
extern int WSInstall P((WSLINK));
extern wsapi_packet WSAnswer P((WSLINK));
extern int WSDoCallPacket P((WSLINK));
#if WSINTERFACE >= 3
extern int WSEvaluate P(( WSLINK, char *));
extern int WSEvaluateString P(( WSLINK, char *));
#else
extern int WSEvaluate P(( WSLINK, charp_ct));
extern int WSEvaluateString P(( WSLINK, charp_ct));
#endif /* WSINTERFACE >= 3 */
WS_END_EXTERN_C

#if WSINTERFACE >= 3
WSMDECL( void, WSDefaultHandler, ( WSLINK, int, int));
#else
WSMDECL( void, WSDefaultHandler, ( WSLINK, unsigned long, unsigned long));
#endif /* WSINTERFACE >= 3 */

#if WSINTERFACE >= 3
WSYDECL( int, WSDefaultYielder, ( WSLINK, WSYieldParameters));
#else
WSYDECL( devyield_result, WSDefaultYielder, ( WSLINK, WSYieldParameters));
#endif /* WSINTERFACE >= 3 */

WS_EXTERN_C
#if WINDOWS_WSTP
extern HWND WSInitializeIcon P(( HINSTANCE hinstCurrent, int nCmdShow));
extern HANDLE WSInstance;
extern HWND WSIconWindow;
#endif
extern int WSAbort, WSDone;
extern long WSSpecialCharacter;
WS_END_EXTERN_C

#endif /* _WSTM_H */





#endif /* _WSTP_H */

