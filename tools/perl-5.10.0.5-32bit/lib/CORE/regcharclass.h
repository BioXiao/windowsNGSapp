/*  -*- buffer-read-only: t -*-
 *
 *    regcharclass.h
 *
 *    Copyright (C) 2007, by Larry Wall and others
 *
 *    You may distribute under the terms of either the GNU General Public
 *    License or the Artistic License, as specified in the README file.
 *
 * !!!!!!!   DO NOT EDIT THIS FILE   !!!!!!!
 * This file is built by Porting/regcharclass.pl.
 * 
 * Any changes made here will be lost!
 *
 */

/*
	LNBREAK: Line Break: \R

	"\x0D\x0A"      # CRLF - Network (Windows) line ending
	0x0A            # LF  | LINE FEED
	0x0B            # VT  | VERTICAL TAB
	0x0C            # FF  | FORM FEED
	0x0D            # CR  | CARRIAGE RETURN
	0x85            # NEL | NEXT LINE
	0x2028          # LINE SEPARATOR
	0x2029          # PARAGRAPH SEPARATOR
*/
/*** GENERATED CODE ***/
#define is_LNBREAK(s,is_utf8)                                               \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                        \
: ( 0x0D == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                     \
: ( is_utf8 ) ?                                                             \
    ( ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
    : 0 )                                                                   \
: ( 0x85 == ((U8*)s)[0] ) )

/*** GENERATED CODE ***/
#define is_LNBREAK_safe(s,e,is_utf8)                                        \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                    \
    : ( 0x0D == ((U8*)s)[0] ) ?                                             \
	( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                 \
    : ( is_utf8 ) ?                                                         \
	( ( 0xC2 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                             \
	: ( 0xE2 == ((U8*)s)[0] ) ?                                         \
	    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
	: 0 )                                                               \
    : ( 0x85 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                    \
    : ( 0x0D == ((U8*)s)[0] ) ?                                             \
	( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                 \
    : ( is_utf8 ) ?                                                         \
	( ( ( 0xC2 == ((U8*)s)[0] ) && ( 0x85 == ((U8*)s)[1] ) ) ? 2 : 0 )  \
    : ( 0x85 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 0) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( !( is_utf8 ) ) ?                                                    \
	( 0x85 == ((U8*)s)[0] )                                             \
    : 0 )                                                                   \
: 0 )

/*** GENERATED CODE ***/
#define is_LNBREAK_utf8(s)                                                  \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                        \
: ( 0x0D == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                     \
: ( 0xC2 == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                     \
: ( 0xE2 == ((U8*)s)[0] ) ?                                                 \
    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
: 0 )

/*** GENERATED CODE ***/
#define is_LNBREAK_utf8_safe(s,e)                                           \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                    \
    : ( 0x0D == ((U8*)s)[0] ) ?                                             \
	( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                 \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
    : 0 )                                                                   \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                    \
    : ( 0x0D == ((U8*)s)[0] ) ?                                             \
	( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                 \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : 0 )                                                                   \
: ((e)-(s) > 0) ?                                                           \
    ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D )                          \
: 0 )

/*** GENERATED CODE ***/
#define is_LNBREAK_latin1(s)                                                \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                        \
: ( 0x0D == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                     \
: ( 0x85 == ((U8*)s)[0] ) )

/*** GENERATED CODE ***/
#define is_LNBREAK_latin1_safe(s,e)                                         \
( ((e)-(s) > 1) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0C ) ? 1                    \
    : ( 0x0D == ((U8*)s)[0] ) ?                                             \
	( ( 0x0A == ((U8*)s)[1] ) ? 2 : 1 )                                 \
    : ( 0x85 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 0) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) || 0x85 == ((U8*)s)[0] )\
: 0 )

/*
	HORIZWS: Horizontal Whitespace: \h \H

	0x09            # HT
	0x20            # SPACE
	0xa0            # NBSP
	0x1680          # OGHAM SPACE MARK
	0x180e          # MONGOLIAN VOWEL SEPARATOR
	0x2000          # EN QUAD
	0x2001          # EM QUAD
	0x2002          # EN SPACE
	0x2003          # EM SPACE
	0x2004          # THREE-PER-EM SPACE
	0x2005          # FOUR-PER-EM SPACE
	0x2006          # SIX-PER-EM SPACE
	0x2007          # FIGURE SPACE
	0x2008          # PUNCTUATION SPACE
	0x2009          # THIN SPACE
	0x200A          # HAIR SPACE
	0x202f          # NARROW NO-BREAK SPACE
	0x205f          # MEDIUM MATHEMATICAL SPACE
	0x3000          # IDEOGRAPHIC SPACE
*/
/*** GENERATED CODE ***/
#define is_HORIZWS(s,is_utf8)                                               \
( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                        \
: ( is_utf8 ) ?                                                             \
    ( ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0xA0 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE1 == ((U8*)s)[0] ) ?                                             \
	( ( 0x9A == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x80 == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: ( 0xA0 == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x8E == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: 0 )                                                               \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x80 == ((U8*)s)[1] ) ?                                         \
	    ( ( ( 0x80 <= ((U8*)s)[2] && ((U8*)s)[2] <= 0x8A ) || 0xAF == ((U8*)s)[2] ) ? 3 : 0 )\
	: ( 0x81 == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x9F == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: 0 )                                                               \
    : ( 0xE3 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0x80 == ((U8*)s)[2] ) ) ? 3 : 0 )  \
    : 0 )                                                                   \
: ( 0xA0 == ((U8*)s)[0] ) )

/*** GENERATED CODE ***/
#define is_HORIZWS_safe(s,e,is_utf8)                                        \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                    \
    : ( is_utf8 ) ?                                                         \
	( ( 0xC2 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0xA0 == ((U8*)s)[1] ) ? 2 : 0 )                             \
	: ( 0xE1 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x9A == ((U8*)s)[1] ) ?                                     \
		( ( 0x80 == ((U8*)s)[2] ) ? 3 : 0 )                         \
	    : ( 0xA0 == ((U8*)s)[1] ) ?                                     \
		( ( 0x8E == ((U8*)s)[2] ) ? 3 : 0 )                         \
	    : 0 )                                                           \
	: ( 0xE2 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x80 == ((U8*)s)[1] ) ?                                     \
		( ( ( 0x80 <= ((U8*)s)[2] && ((U8*)s)[2] <= 0x8A ) || 0xAF == ((U8*)s)[2] ) ? 3 : 0 )\
	    : ( 0x81 == ((U8*)s)[1] ) ?                                     \
		( ( 0x9F == ((U8*)s)[2] ) ? 3 : 0 )                         \
	    : 0 )                                                           \
	: ( 0xE3 == ((U8*)s)[0] ) ?                                         \
	    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0x80 == ((U8*)s)[2] ) ) ? 3 : 0 )\
	: 0 )                                                               \
    : ( 0xA0 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                    \
    : ( is_utf8 ) ?                                                         \
	( ( ( 0xC2 == ((U8*)s)[0] ) && ( 0xA0 == ((U8*)s)[1] ) ) ? 2 : 0 )  \
    : ( 0xA0 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 0) ?                                                           \
    ( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                    \
    : ( !( is_utf8 ) ) ?                                                    \
	( 0xA0 == ((U8*)s)[0] )                                             \
    : 0 )                                                                   \
: 0 )

/*** GENERATED CODE ***/
#define is_HORIZWS_utf8(s)                                                  \
( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                        \
: ( 0xC2 == ((U8*)s)[0] ) ?                                                 \
    ( ( 0xA0 == ((U8*)s)[1] ) ? 2 : 0 )                                     \
: ( 0xE1 == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x9A == ((U8*)s)[1] ) ?                                             \
	( ( 0x80 == ((U8*)s)[2] ) ? 3 : 0 )                                 \
    : ( 0xA0 == ((U8*)s)[1] ) ?                                             \
	( ( 0x8E == ((U8*)s)[2] ) ? 3 : 0 )                                 \
    : 0 )                                                                   \
: ( 0xE2 == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x80 == ((U8*)s)[1] ) ?                                             \
	( ( ( 0x80 <= ((U8*)s)[2] && ((U8*)s)[2] <= 0x8A ) || 0xAF == ((U8*)s)[2] ) ? 3 : 0 )\
    : ( 0x81 == ((U8*)s)[1] ) ?                                             \
	( ( 0x9F == ((U8*)s)[2] ) ? 3 : 0 )                                 \
    : 0 )                                                                   \
: ( 0xE3 == ((U8*)s)[0] ) ?                                                 \
    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0x80 == ((U8*)s)[2] ) ) ? 3 : 0 )      \
: 0 )

/*** GENERATED CODE ***/
#define is_HORIZWS_utf8_safe(s,e)                                           \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                    \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0xA0 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE1 == ((U8*)s)[0] ) ?                                             \
	( ( 0x9A == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x80 == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: ( 0xA0 == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x8E == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: 0 )                                                               \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x80 == ((U8*)s)[1] ) ?                                         \
	    ( ( ( 0x80 <= ((U8*)s)[2] && ((U8*)s)[2] <= 0x8A ) || 0xAF == ((U8*)s)[2] ) ? 3 : 0 )\
	: ( 0x81 == ((U8*)s)[1] ) ?                                         \
	    ( ( 0x9F == ((U8*)s)[2] ) ? 3 : 0 )                             \
	: 0 )                                                               \
    : ( 0xE3 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0x80 == ((U8*)s)[2] ) ) ? 3 : 0 )  \
    : 0 )                                                                   \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] ) ? 1                    \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0xA0 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : 0 )                                                                   \
: ((e)-(s) > 0) ?                                                           \
    ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] )                          \
: 0 )

/*** GENERATED CODE ***/
#define is_HORIZWS_latin1(s)                                                \
( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] || 0xA0 == ((U8*)s)[0] )

/*** GENERATED CODE ***/
#define is_HORIZWS_latin1_safe(s,e)                                         \
( ((e)-(s) > 0) ?                                                           \
    ( 0x09 == ((U8*)s)[0] || 0x20 == ((U8*)s)[0] || 0xA0 == ((U8*)s)[0] )   \
: 0 )

/*** GENERATED CODE ***/
#define is_HORIZWS_cp(cp)                                                   \
( 0x09 == cp || ( 0x09 < cp &&                                              \
( 0x20 == cp || ( 0x20 < cp &&                                              \
( 0xA0 == cp || ( 0xA0 < cp &&                                              \
( 0x1680 == cp || ( 0x1680 < cp &&                                          \
( 0x180E == cp || ( 0x180E < cp &&                                          \
( ( 0x2000 <= cp && cp <= 0x200A ) || ( 0x200A < cp &&                      \
( 0x202F == cp || ( 0x202F < cp &&                                          \
( 0x205F == cp || ( 0x205F < cp &&                                          \
0x3000 == cp ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) )

/*
	VERTWS: Vertical Whitespace: \v \V

	0x0A            # LF
	0x0B            # VT
	0x0C            # FF
	0x0D            # CR
	0x85            # NEL
	0x2028          # LINE SEPARATOR
	0x2029          # PARAGRAPH SEPARATOR
*/
/*** GENERATED CODE ***/
#define is_VERTWS(s,is_utf8)                                                \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                        \
: ( is_utf8 ) ?                                                             \
    ( ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
    : 0 )                                                                   \
: ( 0x85 == ((U8*)s)[0] ) )

/*** GENERATED CODE ***/
#define is_VERTWS_safe(s,e,is_utf8)                                         \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( is_utf8 ) ?                                                         \
	( ( 0xC2 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                             \
	: ( 0xE2 == ((U8*)s)[0] ) ?                                         \
	    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
	: 0 )                                                               \
    : ( 0x85 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( is_utf8 ) ?                                                         \
	( ( ( 0xC2 == ((U8*)s)[0] ) && ( 0x85 == ((U8*)s)[1] ) ) ? 2 : 0 )  \
    : ( 0x85 == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 0) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( !( is_utf8 ) ) ?                                                    \
	( 0x85 == ((U8*)s)[0] )                                             \
    : 0 )                                                                   \
: 0 )

/*** GENERATED CODE ***/
#define is_VERTWS_utf8(s)                                                   \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                        \
: ( 0xC2 == ((U8*)s)[0] ) ?                                                 \
    ( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                     \
: ( 0xE2 == ((U8*)s)[0] ) ?                                                 \
    ( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
: 0 )

/*** GENERATED CODE ***/
#define is_VERTWS_utf8_safe(s,e)                                            \
( ((e)-(s) > 2) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xE2 == ((U8*)s)[0] ) ?                                             \
	( ( ( 0x80 == ((U8*)s)[1] ) && ( 0xA8 == ((U8*)s)[2] || 0xA9 == ((U8*)s)[2] ) ) ? 3 : 0 )\
    : 0 )                                                                   \
: ((e)-(s) > 1) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) ? 1                    \
    : ( 0xC2 == ((U8*)s)[0] ) ?                                             \
	( ( 0x85 == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : 0 )                                                                   \
: ((e)-(s) > 0) ?                                                           \
    ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D )                          \
: 0 )

/*** GENERATED CODE ***/
#define is_VERTWS_latin1(s)                                                 \
( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) || 0x85 == ((U8*)s)[0] )

/*** GENERATED CODE ***/
#define is_VERTWS_latin1_safe(s,e)                                          \
( ((e)-(s) > 0) ?                                                           \
    ( ( 0x0A <= ((U8*)s)[0] && ((U8*)s)[0] <= 0x0D ) || 0x85 == ((U8*)s)[0] )\
: 0 )

/*** GENERATED CODE ***/
#define is_VERTWS_cp(cp)                                                    \
( ( 0x0A <= cp && cp <= 0x0D ) || ( 0x0D < cp &&                            \
( 0x85 == cp || ( 0x85 < cp &&                                              \
( 0x2028 == cp || ( 0x2028 < cp &&                                          \
0x2029 == cp ) ) ) ) ) )

/*
	TRICKYFOLD: Problematic fold case letters.

	0x00DF	# LATIN1 SMALL LETTER SHARP S
	0x0390	# GREEK SMALL LETTER IOTA WITH DIALYTIKA AND TONOS
	0x03B0	# GREEK SMALL LETTER UPSILON WITH DIALYTIKA AND TONOS
*/
/*** GENERATED CODE ***/
#define is_TRICKYFOLD(s,is_utf8)                                            \
( ( is_utf8 ) ?                                                             \
    ( ( 0xC3 == ((U8*)s)[0] ) ?                                             \
	( ( 0x9F == ((U8*)s)[1] ) ? 2 : 0 )                                 \
    : ( 0xCE == ((U8*)s)[0] ) ?                                             \
	( ( 0x90 == ((U8*)s)[1] || 0xB0 == ((U8*)s)[1] ) ? 2 : 0 )          \
    : 0 )                                                                   \
: ( 0xDF == ((U8*)s)[0] ) )

/*** GENERATED CODE ***/
#define is_TRICKYFOLD_safe(s,e,is_utf8)                                     \
( ((e)-(s) > 1) ?                                                           \
    ( ( is_utf8 ) ?                                                         \
	( ( 0xC3 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x9F == ((U8*)s)[1] ) ? 2 : 0 )                             \
	: ( 0xCE == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x90 == ((U8*)s)[1] || 0xB0 == ((U8*)s)[1] ) ? 2 : 0 )      \
	: 0 )                                                               \
    : ( 0xDF == ((U8*)s)[0] ) )                                             \
: ((e)-(s) > 0) ?                                                           \
    ( ( !( is_utf8 ) ) ?                                                    \
	( 0xDF == ((U8*)s)[0] )                                             \
    : 0 )                                                                   \
: 0 )

/*** GENERATED CODE ***/
#define is_TRICKYFOLD_cp(cp)                                                \
( 0xDF == cp || ( 0xDF < cp &&                                              \
( 0x390 == cp || ( 0x390 < cp &&                                            \
0x3B0 == cp ) ) ) )

/*** GENERATED CODE ***/
#define what_TRICKYFOLD(s,is_utf8)                                          \
( ( is_utf8 ) ?                                                             \
    ( ( 0xC3 == ((U8*)s)[0] ) ?                                             \
	( ( 0x9F == ((U8*)s)[1] ) ? 0xDF : 0 )                              \
    : ( 0xCE == ((U8*)s)[0] ) ?                                             \
	( ( 0x90 == ((U8*)s)[1] ) ? 0x390                                   \
	: ( 0xB0 == ((U8*)s)[1] ) ? 0x3B0 : 0 )                             \
    : 0 )                                                                   \
: ( 0xDF == ((U8*)s)[0] ) ? 0xDF : 0 )

/*** GENERATED CODE ***/
#define what_TRICKYFOLD_safe(s,e,is_utf8)                                   \
( ((e)-(s) > 1) ?                                                           \
    ( ( is_utf8 ) ?                                                         \
	( ( 0xC3 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x9F == ((U8*)s)[1] ) ? 0xDF : 0 )                          \
	: ( 0xCE == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x90 == ((U8*)s)[1] ) ? 0x390                               \
	    : ( 0xB0 == ((U8*)s)[1] ) ? 0x3B0 : 0 )                         \
	: 0 )                                                               \
    : ( 0xDF == ((U8*)s)[0] ) ? 0xDF : 0 )                                  \
: ((e)-(s) > 0) ?                                                           \
    ( ( ( !( is_utf8 ) ) && ( 0xDF == ((U8*)s)[0] ) ) ? 0xDF : 0 )          \
: 0 )

/*** GENERATED CODE ***/
#define what_len_TRICKYFOLD(s,is_utf8,len)                                  \
( ( is_utf8 ) ?                                                             \
    ( ( 0xC3 == ((U8*)s)[0] ) ?                                             \
	( ( 0x9F == ((U8*)s)[1] ) ? len=2, 0xDF : 0 )                       \
    : ( 0xCE == ((U8*)s)[0] ) ?                                             \
	( ( 0x90 == ((U8*)s)[1] ) ? len=2, 0x390                            \
	: ( 0xB0 == ((U8*)s)[1] ) ? len=2, 0x3B0 : 0 )                      \
    : 0 )                                                                   \
: ( 0xDF == ((U8*)s)[0] ) ? len=1, 0xDF : 0 )

/*** GENERATED CODE ***/
#define what_len_TRICKYFOLD_safe(s,e,is_utf8,len)                           \
( ((e)-(s) > 1) ?                                                           \
    ( ( is_utf8 ) ?                                                         \
	( ( 0xC3 == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x9F == ((U8*)s)[1] ) ? len=2, 0xDF : 0 )                   \
	: ( 0xCE == ((U8*)s)[0] ) ?                                         \
	    ( ( 0x90 == ((U8*)s)[1] ) ? len=2, 0x390                        \
	    : ( 0xB0 == ((U8*)s)[1] ) ? len=2, 0x3B0 : 0 )                  \
	: 0 )                                                               \
    : ( 0xDF == ((U8*)s)[0] ) ? len=1, 0xDF : 0 )                           \
: ((e)-(s) > 0) ?                                                           \
    ( ( ( !( is_utf8 ) ) && ( 0xDF == ((U8*)s)[0] ) ) ? len=1, 0xDF : 0 )   \
: 0 )

/* ex: set ro: */
