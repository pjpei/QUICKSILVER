#pragma once

#include <stdint.h>

#include "driver/vtx/vtx.h"

void serial_msp_vtx_init();
vtx_update_result_t serial_msp_vtx_update();