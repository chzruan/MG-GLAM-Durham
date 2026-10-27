/* shadow/assert.h: assert() for the Rockstar sources compiled into MG-GLAM.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * Found before the system header (this directory comes first on the include path
 * of the Rockstar sources only). A failed assertion goes through Rockstar's own
 * error path (rockstar_abort_impl, error.cpp), which throws rockstar_error; the
 * adapter catches it per FoF group, so a failed check abandons that output's
 * catalogue instead of killing the simulation with SIGABRT. NDEBUG is honoured. */
#ifdef assert
#undef assert
#endif

#ifdef NDEBUG
#define assert(e) ((void)0)
#else
#ifdef __cplusplus
extern "C"
#endif
void rs_assert_fail(const char *expr, const char *file, int line);
#define assert(e) ((e) ? (void)0 : rs_assert_fail(#e, __FILE__, __LINE__))
#endif
