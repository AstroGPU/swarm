/// Fake memory mapping for portability
//
//

#include <stdlib.h>
void *fakemmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
int fakemunmap(void *addr, size_t length);
int fakemsync(void *addr, size_t length, int flags);


#define PROT_READ       0x1             /* Page can be read.  */
#define PROT_WRITE      0x2             /* Page can be written.  */
#define PROT_EXEC       0x4             /* Page can be executed.  */
#define PROT_NONE       0x0             /* Page can not be accessed.  */
#define MAP_FAILED      ((void *) -1)
#define MAP_SHARED      0x01            /* Share changes.  */
#define MAP_PRIVATE     0x02            /* Changes are private.  */
/* Flags to `msync'.  */
#define MS_ASYNC        1               /* Sync memory asynchronously.  */
#define MS_SYNC         4               /* Synchronous memory sync.  */
#define MS_INVALIDATE   2               /* Invalidate the caches.  */

