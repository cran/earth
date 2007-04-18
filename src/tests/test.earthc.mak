# test.earthc.mak: makefile for test.earthc.main.exe with microsoft visual C 6.0

all: test.earthc.res

R_DIR=\a1\r\work

INCL=-I$(R_DIR)\src\include -I.

# Use PROF_FLAGS  if you want to do profiling using profile, prep and plist
# Note: Don't use -Zi and debug:full flags below if you want to do profiling
PROF_FLAGS=-map -mapinfo:exports -mapinfo:fixups -mapinfo:lines -fixed:no
# PROF_FLAGS=

!IF  "$(CFG)" == "Release"
# Fast version (note: I tried -Ox but it appears no faster than -O2 for this code)
# -02 is for fast code
OUTDIR=Release
CFLAGS=-nologo -DSTANDALONE $(RELEASE_BUILD_CFLAGS) -TP -O2 -W3 -ML $(INCL) -Fp$(OUTDIR)\vc60.PCH -Fo"$(OUTDIR)/" -c
LFLAGS=-nologo $(RELEASE_BUILD_LFLAGS) $(PROF_FLAGS)
# To build Rdll.libs see instructions in gnuwin32\README.packages
LIBS=$(R_DIR)\src\gnuwin32\Rdll.lib $(R_DIR)\bin\Rblas.lib
!ENDIF

!IF  "$(CFG)" == "Debug"
# Debugging version
# -Zi is for a debugging build
# -W3 is warning level 3
# -MLd is for the single threaded static debugging runtime library
# -Gr is for fast function calling (can't use because conflicts with GSL lib)
# No need to define _DEBUG, compiler does it for us if -MLd flag is used
# -D_CRTDBG_MAP_ALLOC maps mallocs/frees to their debug version
OUTDIR=Debug
CFLAGS=-nologo -DSTANDALONE -D_CRTDBG_MAP_ALLOC -TP -Zi -W3 -MLd $(INCL) -Fp$(OUTDIR)\vc60.PCH -Fo"$(OUTDIR)/" -c
LFLAGS=-nologo -debug:full
LIBS=$(R_DIR)\src\gnuwin32\Rdll.lib $(R_DIR)\bin\Rblas.lib
!ENDIF

OBJ=$(OUTDIR)\earth.obj $(OUTDIR)\test.earthc.obj $(OUTDIR)\lib_xerbla.obj

.c{$(OUTDIR)}.obj::
   cl $(CFLAGS) $< 

clean:
	rm -f test.earthc.main.exe $(OUTDIR)/*.obj $(OUTDIR)/*.res $(OUTDIR)/*.pch *.pdb *.dll *.map *.ilk

test.earthc.main.exe: $(OBJ)
	link $(LFLAGS) -out:test.earthc.main.exe $(OBJ) $(LIBS)

test.earthc.res: test.earthc.main.exe test.earthc.ref
	test.earthc.main.exe > $(OUTDIR)\test.earthc.res
!IF  "$(CFG)" == "Release"
	@echo Following diff should give no output
!ENDIF
!IF  "$(CFG)" == "Debug"
	@echo Following diff may give some differences
!ENDIF
	diff test.earthc.ref $(OUTDIR)\test.earthc.res

$(OUTDIR)/earth.obj: ..\earth.c test.earthc.mak
   cl $(CFLAGS) ..\earth.c

$(OUTDIR)/test.earthc.obj: test.earthc.c ..\earth.c test.earthc.mak

$(OUTDIR)/lib_xerbla.obj: test.earthc.mak
