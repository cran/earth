# test.earthc.mak: makefile for test.earthc.main.exe with Microsoft Visual C 6.0
# This builds the executable, runs it, then diffs the results against the reference.

all: test.earthc.out

R_DIR="%ProgramFiles%\r\R-4.0.2"

INCL=-I$(R_DIR)\src\include -I.

# Use PROF_FLAGS  if you want to do profiling using profile, prep and plist
# Note: Don't use -Zi and -debug flags below if you want to do profiling
# PROF_FLAGS=-map -mapinfo:exports -mapinfo:fixups -mapinfo:lines -fixed:no
PROF_FLAGS=

!IF "$(CFG)" != "Release" && "$(CFG)" != "Debug"
!MESSAGE Invalid configuration "$(CFG)" specified. Use "nmake CFG=Debug" or "nmake CFG=Release". Defaulting to CFG=Debug.
CFG=Debug
!ENDIF

!IF  "$(CFG)" == "Release"
# Fast version (note: I tried -Ox but it appears no faster than -O2 for this code)
# -02 is for fast code
OUTDIR=Release
CFLAGS=-nologo -DSTANDALONE $(RELEASE_BUILD_CFLAGS) -TP -O2 -W3 -MT $(INCL) -Fp$(OUTDIR)\vc60.PCH -Fo"$(OUTDIR)/" -c
LFLAGS=-nologo $(RELEASE_BUILD_LFLAGS) $(PROF_FLAGS)
# To build R.lib, see for example https://www.asawicki.info/news_1420_generating_lib_file_for_dll_library
# See also c:\bin\ddl2lib.bat
LIBS=R.lib Rblas.lib
!ENDIF

!IF  "$(CFG)" == "Debug"
# Debugging version
# -Tp says treat the file as a C++ file (needed for C99 source files)
# -Zi is for a debugging build
# -W3 is warning level 3
# -MTd is for the multi threaded static debugging runtime library
# -Gr is for fast function calling (can't use because conflicts with GSL lib)
# No need to define _DEBUG, compiler does it for us if -MTd flag is used
OUTDIR=Debug
# TODO We use -MDd instead of -MTd here.  It seems to work just as well.
#      Using -Mtd causes linker error: LIBCMTD.lib(sprintf.obj) : error LNK2005: _sprintf already defined in R.lib(R.dll)
# CFLAGS=-nologo -DSTANDALONE -TP -Zi -W3 -MTd $(INCL) -Fp$(OUTDIR)\vc60.PCH -Fo"$(OUTDIR)/" -c
CFLAGS=-nologo -DSTANDALONE -TP -Zi -W3 -MDd $(INCL) -Fp$(OUTDIR)\vc60.PCH -Fo"$(OUTDIR)/" -c
LFLAGS=-nologo -debug
# To build R.libs see instructions in gnuwin32\README.packages
LIBS=R.lib Rblas.lib 
!ENDIF

OBJ=$(OUTDIR)\earth.obj $(OUTDIR)\test.earthc.obj

.c{$(OUTDIR)}.obj::
   cl $(CFLAGS) $< 

clean:
	rm -f test.earthc.main.exe $(OUTDIR)/*.obj $(OUTDIR)/*.out $(OUTDIR)/*.pch *.pdb *.dll *.map *.ilk

test.earthc.main.exe: $(OBJ)
	@rem Added path below to disambiguate from rtools which (may 2020 R version 4.0.0)
	"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\BIN\link.exe" $(LFLAGS) -out:test.earthc.main.exe $(OBJ) $(LIBS)

# we use diff -w below so \r\n is treated the same as \n
test.earthc.out: test.earthc.main.exe test.earthc.out.save
	@rem following is so we can check the compiler version because it can affect the model
	-cl 1> NUL 2> $(OUTDIR)\test.earthc.out
	test.earthc.main.exe >> $(OUTDIR)\test.earthc.out
!IF  "$(CFG)" == "Debug"
	@rem @echo === Following diff may give some differences ===
!ENDIF
	diff -w $(OUTDIR)\test.earthc.out test.earthc.out.save 

$(OUTDIR)/earth.obj: ..\..\src\earth.c test.earthc.mak
   cl $(CFLAGS) ..\..\src\earth.c

$(OUTDIR)/test.earthc.obj: test.earthc.c ..\..\src\earth.c test.earthc.mak
