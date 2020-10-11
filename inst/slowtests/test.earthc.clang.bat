@rem test.earthc.clang.bat: test the standalone earth.c with main()
@rem
@rem This tests the earth C code.  It does this: builds test.earthc.exe
@rem (under clang), runs it, and compares results to test.earthc.out.save.
@rem
@rem You will need to tweak this file and test.earthc.gcc.mak for your directories.
@rem
@rem You need to make R.lib first -- see instructions in gnuwin32/README.packages.

@echo test.earthc.clang.bat
@mks.cp "C:\bin\R400devdll\i386\R.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create R.lib and Rblas.lib beforehand
@mks.cp "C:\bin\R400devdll\i386\R.lib" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rblas.lib" .
                                @if %errorlevel% neq 0 goto error

@set CLANGEXE="c:\Program Files (x86)/LLVM/bin/clang.exe"

@rem Flags same as gcc where possible.
@rem -U_MSC_VER is needed because some clang executables define this inherently https://stackoverflow.com/questions/38499462/how-to-tell-clang-to-stop-pretending-to-be-other-compilers
@rem Some of these warning suppressions are necessary because we use R and BLAS routines.
@rem Oct 2020: tested with clang version 6.0.0 (tags/RELEASE_600/final)

%CLANGEXE% -DSTANDALONE -Wall --pedantic -Wextra -Weverything -O3 -std=gnu99^
 -U_MSC_VER^
 -Wno-strict-prototypes -Wno-reserved-id-macro -Wno-cast-qual -Wno-unknown-pragmas^
 -Wno-float-equal -Wno-format-nonliteral -Wno-padded -Wno-sign-conversion -Wno-undef^
 -Wno-shadow -Wno-deprecated-declarations -Wno-implicit-function-declaration^
 -Wno-missing-noreturn -Wno-missing-prototypes -Wno-unused-parameter^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c test.earthc.c^
 R.lib Rblas.lib -o earthc-clang.exe
                                @if %errorlevel% neq 0 goto error
earthc-clang.exe > test.earthc-clang.out
                                @rem no errorlevel test, diff will do check for discrepancies
                                @rem @if %errorlevel% neq 0 goto error
mks.diff test.earthc-clang.out test.earthc.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll R.lib Rblas.lib earthc-clang.* test.earthc-clang.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
