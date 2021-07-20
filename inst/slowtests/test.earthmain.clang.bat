@rem test.earthmain.clang.bat: test the standalone earth.c with main()
@rem
@rem The gcc, Microsoft, and clang compiler batch files all test
@rem against the same reference file "test.earthmain.out.save"
@rem
@rem Stephen Milborrow Dec 2014 Shrewsbury.  Updated Petaluma May 2020.

@echo test.earthmain.clang.bat

@rem The following is a basic check that you have Visual Studio 2019 for 32 bit targets
@which cl | egrep -i "Visual.Studio.2019.Community.VC.Tools.MSVC.*.bin.Hostx.*x86.cl" >NUL && goto donesetpath
@echo Environment is not VC16 (Visual Studio 2019) 32 bit -- please invoke vc16-32.bat
@exit /B 1
:donesetpath

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
%CLANGEXE% --version

%CLANGEXE% -DSTANDALONE -DMAIN -Wall --pedantic -Wextra -Weverything -O3 -std=gnu99^
 -U_MSC_VER^
 -Wno-strict-prototypes -Wno-reserved-id-macro -Wno-cast-qual -Wno-unknown-pragmas^
 -Wno-float-equal -Wno-format-nonliteral -Wno-padded -Wno-sign-conversion -Wno-undef^
 -Wno-shadow -Wno-deprecated-declarations -Wno-implicit-function-declaration^
 -Wno-missing-noreturn^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 R.lib Rblas.lib -o earthmain-clang.exe
                                @if %errorlevel% neq 0 goto error
earthmain-clang.exe > test.earthmain-clang.out
                                @rem no errorlevel test, diff will do check for discrepancies
                                @rem @if %errorlevel% neq 0 goto error
mks.diff test.earthmain-clang.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll R.lib Rblas.lib earthmain-clang.* test.earthmain-clang.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
