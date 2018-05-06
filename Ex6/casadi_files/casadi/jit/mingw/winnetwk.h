/**
 * This file is part of the mingw-w64 runtime package.
 * No warranty is given; refer to the file DISCLAIMER within this package.
 */
#ifndef _WINNETWK_
#define _WINNETWK_

#include <winapifamily.h>
#if WINAPI_FAMILY_PARTITION(WINAPI_PARTITION_DESKTOP)

#include <_mingw_unicode.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <wnnc.h>

#define RESOURCE_CONNECTED 0x00000001
#define RESOURCE_GLOBALNET 0x00000002
#define RESOURCE_REMEMBERED 0x00000003
#define RESOURCE_RECENT 0x00000004
#define RESOURCE_CONTEXT 0x00000005

#define RESOURCETYPE_ANY 0x00000000
#define RESOURCETYPE_DISK 0x00000001
#define RESOURCETYPE_PRINT 0x00000002
#define RESOURCETYPE_RESERVED 0x00000008

#define RESOURCETYPE_UNKNOWN 0xFFFFFFFF

#define RESOURCEUSAGE_CONNECTABLE 0x00000001
#define RESOURCEUSAGE_CONTAINER 0x00000002
#define RESOURCEUSAGE_NOLOCALDEVICE 0x00000004
#define RESOURCEUSAGE_SIBLING 0x00000008
#define RESOURCEUSAGE_ATTACHED 0x00000010
#define RESOURCEUSAGE_ALL (RESOURCEUSAGE_CONNECTABLE | RESOURCEUSAGE_CONTAINER | RESOURCEUSAGE_ATTACHED)
#define RESOURCEUSAGE_RESERVED 0x80000000

#define RESOURCEDISPLAYTYPE_GENERIC 0x00000000
#define RESOURCEDISPLAYTYPE_DOMAIN 0x00000001
#define RESOURCEDISPLAYTYPE_SERVER 0x00000002
#define RESOURCEDISPLAYTYPE_SHARE 0x00000003
#define RESOURCEDISPLAYTYPE_FILE 0x00000004
#define RESOURCEDISPLAYTYPE_GROUP 0x00000005
#define RESOURCEDISPLAYTYPE_NETWORK 0x00000006
#define RESOURCEDISPLAYTYPE_ROOT 0x00000007
#define RESOURCEDISPLAYTYPE_SHAREADMIN 0x00000008
#define RESOURCEDISPLAYTYPE_DIRECTORY 0x00000009
#define RESOURCEDISPLAYTYPE_TREE 0x0000000a
#define RESOURCEDISPLAYTYPE_NDSCONTAINER 0x0000000b

  typedef struct _NETRESOURCEA {
    DWORD dwScope;
    DWORD dwType;
    DWORD dwDisplayType;
    DWORD dwUsage;
    LPSTR lpLocalName;
    LPSTR lpRemoteName;
    LPSTR lpComment;
    LPSTR lpProvider;
  } NETRESOURCEA,*LPNETRESOURCEA;
  typedef struct _NETRESOURCEW {
    DWORD dwScope;
    DWORD dwType;
    DWORD dwDisplayType;
    DWORD dwUsage;
    LPWSTR lpLocalName;
    LPWSTR lpRemoteName;
    LPWSTR lpComment;
    LPWSTR lpProvider;
  } NETRESOURCEW,*LPNETRESOURCEW;

  __MINGW_TYPEDEF_AW(NETRESOURCE)
  __MINGW_TYPEDEF_AW(LPNETRESOURCE)

#define NETPROPERTY_PERSISTENT 1

#define CONNECT_UPDATE_PROFILE 0x00000001
#define CONNECT_UPDATE_RECENT 0x00000002
#define CONNECT_TEMPORARY 0x00000004
#define CONNECT_INTERACTIVE 0x00000008
#define CONNECT_PROMPT 0x00000010
#define CONNECT_NEED_DRIVE 0x00000020
#define CONNECT_REFCOUNT 0x00000040
#define CONNECT_REDIRECT 0x00000080
#define CONNECT_LOCALDRIVE 0x00000100
#define CONNECT_CURRENT_MEDIA 0x00000200
#define CONNECT_DEFERRED 0x00000400
#define CONNECT_RESERVED 0xFF000000
#define CONNECT_COMMANDLINE 0x00000800
#define CONNECT_CMD_SAVECRED 0x00001000
#if WINVER >= 0x0600
#define CONNECT_CRED_RESET 0x00002000
#endif

#define WNetAddConnection __MINGW_NAME_AW(WNetAddConnection)
#define WNetAddConnection2 __MINGW_NAME_AW(WNetAddConnection2)
#define WNetAddConnection3 __MINGW_NAME_AW(WNetAddConnection3)
#define WNetCancelConnection __MINGW_NAME_AW(WNetCancelConnection)
#define WNetCancelConnection2 __MINGW_NAME_AW(WNetCancelConnection2)
#define WNetGetConnection __MINGW_NAME_AW(WNetGetConnection)
#define WNetRestoreConnection __MINGW_NAME_AW(WNetRestoreConnection)
#define WNetUseConnection __MINGW_NAME_AW(WNetUseConnection)

  DWORD WINAPI WNetAddConnectionA(LPCSTR lpRemoteName,LPCSTR lpPassword,LPCSTR lpLocalName);
  DWORD WINAPI WNetAddConnectionW(LPCWSTR lpRemoteName,LPCWSTR lpPassword,LPCWSTR lpLocalName);
  DWORD WINAPI WNetAddConnection2A(LPNETRESOURCEA lpNetResource,LPCSTR lpPassword,LPCSTR lpUserName,DWORD dwFlags);
  DWORD WINAPI WNetAddConnection2W(LPNETRESOURCEW lpNetResource,LPCWSTR lpPassword,LPCWSTR lpUserName,DWORD dwFlags);
  DWORD WINAPI WNetAddConnection3A(HWND hwndOwner,LPNETRESOURCEA lpNetResource,LPCSTR lpPassword,LPCSTR lpUserName,DWORD dwFlags);
  DWORD WINAPI WNetAddConnection3W(HWND hwndOwner,LPNETRESOURCEW lpNetResource,LPCWSTR lpPassword,LPCWSTR lpUserName,DWORD dwFlags);
  DWORD WINAPI WNetCancelConnectionA(LPCSTR lpName,WINBOOL fForce);
  DWORD WINAPI WNetCancelConnectionW(LPCWSTR lpName,WINBOOL fForce);
  DWORD WINAPI WNetCancelConnection2A(LPCSTR lpName,DWORD dwFlags,WINBOOL fForce);
  DWORD WINAPI WNetCancelConnection2W(LPCWSTR lpName,DWORD dwFlags,WINBOOL fForce);
  DWORD WINAPI WNetGetConnectionA(LPCSTR lpLocalName,LPSTR lpRemoteName,LPDWORD lpnLength);
  DWORD WINAPI WNetGetConnectionW(LPCWSTR lpLocalName,LPWSTR lpRemoteName,LPDWORD lpnLength);
  DWORD WINAPI WNetRestoreConnectionA(HWND hwndParent,LPCSTR lpDevice);
  DWORD WINAPI WNetUseConnectionA(HWND hwndOwner,LPNETRESOURCEA lpNetResource,LPCSTR lpPassword,LPCSTR lpUserID,DWORD dwFlags,LPSTR lpAccessName,LPDWORD lpBufferSize,LPDWORD lpResult);
  DWORD WINAPI WNetUseConnectionW(HWND hwndOwner,LPNETRESOURCEW lpNetResource,LPCWSTR lpPassword,LPCWSTR lpUserID,DWORD dwFlags,LPWSTR lpAccessName,LPDWORD lpBufferSize,LPDWORD lpResult);
  DWORD WINAPI WNetConnectionDialog(HWND hwnd,DWORD dwType);
  DWORD WINAPI WNetDisconnectDialog(HWND hwnd,DWORD dwType);
#if (_WIN32_WINNT >= _WIN32_WINNT_LONGHORN)
  DWORD WINAPI WNetRestoreSingleConnectionW(HWND hwndParent, LPCWSTR lpDevice, BOOL fUseUI);
#else
  DWORD WINAPI WNetRestoreConnectionW(HWND hwndParent,LPCWSTR lpDevice);
#endif

  typedef struct _CONNECTDLGSTRUCTA {
    DWORD cbStructure;
    HWND hwndOwner;
    LPNETRESOURCEA lpConnRes;
    DWORD dwFlags;
    DWORD dwDevNum;
  } CONNECTDLGSTRUCTA,*LPCONNECTDLGSTRUCTA;

  typedef struct _CONNECTDLGSTRUCTW {
    DWORD cbStructure;
    HWND hwndOwner;
    LPNETRESOURCEW lpConnRes;
    DWORD dwFlags;
    DWORD dwDevNum;
  } CONNECTDLGSTRUCTW,*LPCONNECTDLGSTRUCTW;

  __MINGW_TYPEDEF_AW(CONNECTDLGSTRUCT)
  __MINGW_TYPEDEF_AW(LPCONNECTDLGSTRUCT)

#define CONNDLG_RO_PATH 0x00000001
#define CONNDLG_CONN_POINT 0x00000002
#define CONNDLG_USE_MRU 0x00000004
#define CONNDLG_HIDE_BOX 0x00000008

#define CONNDLG_PERSIST 0x00000010
#define CONNDLG_NOT_PERSIST 0x00000020

#define WNetConnectionDialog1 __MINGW_NAME_AW(WNetConnectionDialog1)

  DWORD WINAPI WNetConnectionDialog1A(LPCONNECTDLGSTRUCTA lpConnDlgStruct);
  DWORD WINAPI WNetConnectionDialog1W(LPCONNECTDLGSTRUCTW lpConnDlgStruct);

  typedef struct _DISCDLGSTRUCTA {
    DWORD cbStructure;
    HWND hwndOwner;
    LPSTR lpLocalName;
    LPSTR lpRemoteName;
    DWORD dwFlags;
  } DISCDLGSTRUCTA,*LPDISCDLGSTRUCTA;

  typedef struct _DISCDLGSTRUCTW {
    DWORD cbStructure;
    HWND hwndOwner;
    LPWSTR lpLocalName;
    LPWSTR lpRemoteName;
    DWORD dwFlags;
  } DISCDLGSTRUCTW,*LPDISCDLGSTRUCTW;

  __MINGW_TYPEDEF_AW(DISCDLGSTRUCT)
  __MINGW_TYPEDEF_AW(LPDISCDLGSTRUCT)

#define DISC_UPDATE_PROFILE 0x00000001
#define DISC_NO_FORCE 0x00000040

#define WNetDisconnectDialog1 __MINGW_NAME_AW(WNetDisconnectDialog1)
#define WNetOpenEnum __MINGW_NAME_AW(WNetOpenEnum)
#define WNetEnumResource __MINGW_NAME_AW(WNetEnumResource)
#define WNetGetResourceParent __MINGW_NAME_AW(WNetGetResourceParent)
#define WNetGetResourceInformation __MINGW_NAME_AW(WNetGetResourceInformation)

  DWORD WINAPI WNetDisconnectDialog1A(LPDISCDLGSTRUCTA lpConnDlgStruct);
  DWORD WINAPI WNetDisconnectDialog1W(LPDISCDLGSTRUCTW lpConnDlgStruct);
  DWORD WINAPI WNetOpenEnumA(DWORD dwScope,DWORD dwType,DWORD dwUsage,LPNETRESOURCEA lpNetResource,LPHANDLE lphEnum);
  DWORD WINAPI WNetOpenEnumW(DWORD dwScope,DWORD dwType,DWORD dwUsage,LPNETRESOURCEW lpNetResource,LPHANDLE lphEnum);
  DWORD WINAPI WNetEnumResourceA(HANDLE hEnum,LPDWORD lpcCount,LPVOID lpBuffer,LPDWORD lpBufferSize);
  DWORD WINAPI WNetEnumResourceW(HANDLE hEnum,LPDWORD lpcCount,LPVOID lpBuffer,LPDWORD lpBufferSize);
  DWORD WINAPI WNetCloseEnum(HANDLE hEnum);
  DWORD WINAPI WNetGetResourceParentA(LPNETRESOURCEA lpNetResource,LPVOID lpBuffer,LPDWORD lpcbBuffer);
  DWORD WINAPI WNetGetResourceParentW(LPNETRESOURCEW lpNetResource,LPVOID lpBuffer,LPDWORD lpcbBuffer);
  DWORD WINAPI WNetGetResourceInformationA(LPNETRESOURCEA lpNetResource,LPVOID lpBuffer,LPDWORD lpcbBuffer,LPSTR *lplpSystem);
  DWORD WINAPI WNetGetResourceInformationW(LPNETRESOURCEW lpNetResource,LPVOID lpBuffer,LPDWORD lpcbBuffer,LPWSTR *lplpSystem);

#define UNIVERSAL_NAME_INFO_LEVEL 0x00000001
#define REMOTE_NAME_INFO_LEVEL 0x00000002

  typedef struct _UNIVERSAL_NAME_INFOA {
    LPSTR lpUniversalName;
  } UNIVERSAL_NAME_INFOA,*LPUNIVERSAL_NAME_INFOA;

  typedef struct _UNIVERSAL_NAME_INFOW {
    LPWSTR lpUniversalName;
  } UNIVERSAL_NAME_INFOW,*LPUNIVERSAL_NAME_INFOW;

  __MINGW_TYPEDEF_AW(UNIVERSAL_NAME_INFO)
  __MINGW_TYPEDEF_AW(LPUNIVERSAL_NAME_INFO)

  typedef struct _REMOTE_NAME_INFOA {
    LPSTR lpUniversalName;
    LPSTR lpConnectionName;
    LPSTR lpRemainingPath;
  } REMOTE_NAME_INFOA,*LPREMOTE_NAME_INFOA;

  typedef struct _REMOTE_NAME_INFOW {
    LPWSTR lpUniversalName;
    LPWSTR lpConnectionName;
    LPWSTR lpRemainingPath;
  } REMOTE_NAME_INFOW,*LPREMOTE_NAME_INFOW;

  __MINGW_TYPEDEF_AW(REMOTE_NAME_INFO)
  __MINGW_TYPEDEF_AW(LPREMOTE_NAME_INFO)

#define WNetGetUniversalName __MINGW_NAME_AW(WNetGetUniversalName)
#define WNetGetUser __MINGW_NAME_AW(WNetGetUser)
#define WNetGetProviderName __MINGW_NAME_AW(WNetGetProviderName)

  DWORD WINAPI WNetGetUniversalNameA(LPCSTR lpLocalPath,DWORD dwInfoLevel,LPVOID lpBuffer,LPDWORD lpBufferSize);
  DWORD WINAPI WNetGetUniversalNameW(LPCWSTR lpLocalPath,DWORD dwInfoLevel,LPVOID lpBuffer,LPDWORD lpBufferSize);
  DWORD WINAPI WNetGetUserA(LPCSTR lpName,LPSTR lpUserName,LPDWORD lpnLength);
  DWORD WINAPI WNetGetUserW(LPCWSTR lpName,LPWSTR lpUserName,LPDWORD lpnLength);

#define WNFMT_MULTILINE 0x01
#define WNFMT_ABBREVIATED 0x02
#define WNFMT_INENUM 0x10
#define WNFMT_CONNECTION 0x20

  DWORD WINAPI WNetGetProviderNameA(DWORD dwNetType,LPSTR lpProviderName,LPDWORD lpBufferSize);
  DWORD WINAPI WNetGetProviderNameW(DWORD dwNetType,LPWSTR lpProviderName,LPDWORD lpBufferSize);

  typedef struct _NETINFOSTRUCT {
    DWORD cbStructure;
    DWORD dwProviderVersion;
    DWORD dwStatus;
    DWORD dwCharacteristics;
    ULONG_PTR dwHandle;
    WORD wNetType;
    DWORD dwPrinters;
    DWORD dwDrives;
  } NETINFOSTRUCT,*LPNETINFOSTRUCT;

#define NETINFO_DLL16 0x00000001
#define NETINFO_DISKRED 0x00000004
#define NETINFO_PRINTERRED 0x00000008

#define WNetGetNetworkInformation __MINGW_NAME_AW(WNetGetNetworkInformation)

  DWORD WINAPI WNetGetNetworkInformationA(LPCSTR lpProvider,LPNETINFOSTRUCT lpNetInfoStruct);
  DWORD WINAPI WNetGetNetworkInformationW(LPCWSTR lpProvider,LPNETINFOSTRUCT lpNetInfoStruct);

  typedef UINT (WINAPI *PFNGETPROFILEPATHA) (LPCSTR pszUsername,LPSTR pszBuffer,UINT cbBuffer);
  typedef UINT (WINAPI *PFNGETPROFILEPATHW) (LPCWSTR pszUsername,LPWSTR pszBuffer,UINT cbBuffer);

#define PFNGETPROFILEPATH __MINGW_NAME_AW(PFNGETPROFILEPATH)

  typedef UINT (WINAPI *PFNRECONCILEPROFILEA) (LPCSTR pszCentralFile,LPCSTR pszLocalFile,DWORD dwFlags);
  typedef UINT (WINAPI *PFNRECONCILEPROFILEW) (LPCWSTR pszCentralFile,LPCWSTR pszLocalFile,DWORD dwFlags);

#define PFNRECONCILEPROFILE __MINGW_NAME_AW(PFNRECONCILEPROFILE)

#define RP_LOGON 0x01
#define RP_INIFILE 0x02

  typedef WINBOOL (WINAPI *PFNPROCESSPOLICIESA) (HWND hwnd,LPCSTR pszPath,LPCSTR pszUsername,LPCSTR pszComputerName,DWORD dwFlags);
  typedef WINBOOL (WINAPI *PFNPROCESSPOLICIESW) (HWND hwnd,LPCWSTR pszPath,LPCWSTR pszUsername,LPCWSTR pszComputerName,DWORD dwFlags);

#define PFNPROCESSPOLICIES __MINGW_NAME_AW(PFNPROCESSPOLICIES)

#define PP_DISPLAYERRORS 0x01

#define WNetGetLastError __MINGW_NAME_AW(WNetGetLastError)

  DWORD WINAPI WNetGetLastErrorA(LPDWORD lpError,LPSTR lpErrorBuf,DWORD nErrorBufSize,LPSTR lpNameBuf,DWORD nNameBufSize);
  DWORD WINAPI WNetGetLastErrorW(LPDWORD lpError,LPWSTR lpErrorBuf,DWORD nErrorBufSize,LPWSTR lpNameBuf,DWORD nNameBufSize);

#define WN_SUCCESS NO_ERROR
#define WN_NO_ERROR NO_ERROR
#define WN_NOT_SUPPORTED ERROR_NOT_SUPPORTED
#define WN_CANCEL ERROR_CANCELLED
#define WN_RETRY ERROR_RETRY
#define WN_NET_ERROR ERROR_UNEXP_NET_ERR
#define WN_MORE_DATA ERROR_MORE_DATA
#define WN_BAD_POINTER ERROR_INVALID_ADDRESS
#define WN_BAD_VALUE ERROR_INVALID_PARAMETER
#define WN_BAD_USER ERROR_BAD_USERNAME
#define WN_BAD_PASSWORD ERROR_INVALID_PASSWORD
#define WN_ACCESS_DENIED ERROR_ACCESS_DENIED
#define WN_FUNCTION_BUSY ERROR_BUSY
#define WN_WINDOWS_ERROR ERROR_UNEXP_NET_ERR
#define WN_OUT_OF_MEMORY ERROR_NOT_ENOUGH_MEMORY
#define WN_NO_NETWORK ERROR_NO_NETWORK
#define WN_EXTENDED_ERROR ERROR_EXTENDED_ERROR
#define WN_BAD_LEVEL ERROR_INVALID_LEVEL
#define WN_BAD_HANDLE ERROR_INVALID_HANDLE
#define WN_NOT_INITIALIZING ERROR_ALREADY_INITIALIZED
#define WN_NO_MORE_DEVICES ERROR_NO_MORE_DEVICES
#define WN_NOT_CONNECTED ERROR_NOT_CONNECTED
#define WN_OPEN_FILES ERROR_OPEN_FILES
#define WN_DEVICE_IN_USE ERROR_DEVICE_IN_USE
#define WN_BAD_NETNAME ERROR_BAD_NET_NAME
#define WN_BAD_LOCALNAME ERROR_BAD_DEVICE
#define WN_ALREADY_CONNECTED ERROR_ALREADY_ASSIGNED
#define WN_DEVICE_ERROR ERROR_GEN_FAILURE
#define WN_CONNECTION_CLOSED ERROR_CONNECTION_UNAVAIL
#define WN_NO_NET_OR_BAD_PATH ERROR_NO_NET_OR_BAD_PATH
#define WN_BAD_PROVIDER ERROR_BAD_PROVIDER
#define WN_CANNOT_OPEN_PROFILE ERROR_CANNOT_OPEN_PROFILE
#define WN_BAD_PROFILE ERROR_BAD_PROFILE
#define WN_BAD_DEV_TYPE ERROR_BAD_DEV_TYPE
#define WN_DEVICE_ALREADY_REMEMBERED ERROR_DEVICE_ALREADY_REMEMBERED
#define WN_CONNECTED_OTHER_PASSWORD ERROR_CONNECTED_OTHER_PASSWORD
#define WN_CONNECTED_OTHER_PASSWORD_DEFAULT ERROR_CONNECTED_OTHER_PASSWORD_DEFAULT
#define WN_NO_MORE_ENTRIES ERROR_NO_MORE_ITEMS
#define WN_NOT_CONTAINER ERROR_NOT_CONTAINER
#define WN_NOT_AUTHENTICATED ERROR_NOT_AUTHENTICATED
#define WN_NOT_LOGGED_ON ERROR_NOT_LOGGED_ON
#define WN_NOT_VALIDATED ERROR_NO_LOGON_SERVERS

  typedef struct _NETCONNECTINFOSTRUCT {
    DWORD cbStructure;
    DWORD dwFlags;
    DWORD dwSpeed;
    DWORD dwDelay;
    DWORD dwOptDataSize;
  } NETCONNECTINFOSTRUCT,*LPNETCONNECTINFOSTRUCT;

#define WNCON_FORNETCARD 0x00000001
#define WNCON_NOTROUTED 0x00000002
#define WNCON_SLOWLINK 0x00000004
#define WNCON_DYNAMIC 0x00000008

#define MultinetGetConnectionPerformance __MINGW_NAME_AW(MultinetGetConnectionPerformance)

  DWORD WINAPI MultinetGetConnectionPerformanceA(LPNETRESOURCEA lpNetResource,LPNETCONNECTINFOSTRUCT lpNetConnectInfoStruct);
  DWORD WINAPI MultinetGetConnectionPerformanceW(LPNETRESOURCEW lpNetResource,LPNETCONNECTINFOSTRUCT lpNetConnectInfoStruct);

#ifdef __cplusplus
}
#endif

#endif /* WINAPI_PARTITION_DESKTOP.  */

#endif
