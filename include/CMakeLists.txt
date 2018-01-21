CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/fteik_config.h)
SET(headers fteik_config.h fteik_fortran.h fteik.h fteik_h5io.h fteik_io.h
            fteik_os.h fteik_sorting.h fteik_struct.h fteik_xdmf.h)
FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
INSTALL(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
