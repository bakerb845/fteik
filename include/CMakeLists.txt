configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/fteik_config.h)
set(headers fteik_config.h fteik.h fteik_fortran.h fteik_sorting.h)
file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
install(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
