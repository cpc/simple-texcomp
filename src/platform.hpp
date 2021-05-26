#ifndef PLATFORM_H
#define PLATFORM_H

/*
 * Platform-specific definitions
 */

#if __ANDROID__

#include <android/log.h>
#define LOGI(...) \
    ((void)__android_log_print(ANDROID_LOG_INFO, "native-activity", __VA_ARGS__))
#define LOGW(...) \
    ((void)__android_log_print(ANDROID_LOG_WARN, "native-activity", __VA_ARGS__))
#define LOGE(...) \
    ((void)__android_log_print(ANDROID_LOG_ERROR, "native-activity", __VA_ARGS__))

#else // __ANDROID__

#include <cstdio>
#define LOGI(...) printf(__VA_ARGS__)
#define LOGW(...) printf("WARNING: " __VA_ARGS__)
#define LOGE(...) printf("ERROR: " __VA_ARGS__)

#endif // __ANDROID__

#endif // PLATFORM_H
