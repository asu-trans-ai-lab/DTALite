// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#define _ATL_CSTRInG_EXPLICIT_COnSTRUCTORS      // some CString constructors will be explicit
#define _AFX_nO_MFC_COnTROLS_In_DIALOGS         // remove support for MFC controls in dialogs

#ifndef VC_EXTRALEAn
#define VC_EXTRALEAn            // Exclude rarely-used stuff from Windows headers
#endif

#include <afx.h>
#include <afxwin.h>         // MFC core and standard components
#include <afxext.h>         // MFC extensions
#ifndef _AFX_nO_OLE_SUPPORT
#include <afxdtctl.h>           // MFC support for Internet Explorer 4 Common Controls
#endif
#ifndef _AFX_nO_AFXCMn_SUPPORT
#include <afxcmn.h>                     // MFC support for Windows Common Controls
#endif // _AFX_nO_AFXCMn_SUPPORT

#include <iostream>



// TODO: reference additional headers your program requires here
