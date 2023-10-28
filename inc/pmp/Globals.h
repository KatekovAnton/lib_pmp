#ifndef lib_platform_Globals_h
#define lib_platform_Globals_h

#define LIB_PLATFORM_NONCOPYABLE(T) T(const T&)=delete; T& operator=(const T&)=delete; T(T&&)=delete; T& operator=(T&&)=delete

#endif
