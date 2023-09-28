/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_PROPERTY(density, cell, thread);
extern DEFINE_PROPERTY(viscosity, cell, thread);
extern DEFINE_PROPERTY(speed_sound, cell, thread);
extern DEFINE_EXECUTE_AT_END(eval_barotropic);
__declspec(dllexport) UDF_Data udf_data[] = {
{"density", (void(*)())density, UDF_TYPE_PROPERTY},
{"viscosity", (void(*)())viscosity, UDF_TYPE_PROPERTY},
{"speed_sound", (void(*)())speed_sound, UDF_TYPE_PROPERTY},
{"eval_barotropic", (void(*)())eval_barotropic, UDF_TYPE_EXECUTE_AT_END},
};
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
