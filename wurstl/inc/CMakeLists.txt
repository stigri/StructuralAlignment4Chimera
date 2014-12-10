file(GLOB ALL_FILES *)
foreach(full_name ${ALL_FILES})
  get_filename_component(fname ${full_name} NAME)
  if(${fname} MATCHES ".h")
    set(INC_HEADERS
    ${INC_HEADERS}
    ${fname}
    )
  endif()
endforeach()


