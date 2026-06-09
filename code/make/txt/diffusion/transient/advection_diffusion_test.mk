.PHONY: txt_advection_diffusion_test\
        txt_advection_diffusion_test_clean\
        txt_advection_diffusion_test_distclean

build/$(MODE)/txt/diffusion/transient/advection_diffusion_test:\
  | build/$(MODE)/txt/diffusion/transient
	mkdir -p $@

_txt_diffusion_transient_advection_diffusion_test :=\
  build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/mixed_weak_cochain_brick_3d_2_forman_input.txt\
  build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/primal_weak_cochain_brick_3d_2_forman_input.txt\
  build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/primal_weak_cochain_brick_3d_2_forman_trapezoidal_0p01_200_potential.txt

build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/mixed_weak_cochain_brick_3d_2_forman_input.txt:\
  build/$(MODE)/bin/diffusion_transient_discrete_mixed_weak_from_continuous$(.EXE)\
  build/$(MODE)/txt/mesh/brick_3d_2/forman.txt\
  build/$(MODE)/txt/mesh/brick_3d_2/forman_hodge_corrected.txt\
  build/$(MODE)/txt/mesh/brick_3d_2/forman_vol.txt\
  build/$(MODE)/obj/plugins/diffusion_transient_advection_diffusion_test_fluid_flow$(.OBJ)\
  | build/$(MODE)/txt/diffusion/transient/advection_diffusion_test\
    build/$(MODE)/lib/plugins/libdiffusion$(.DLL)
	$(INTERPRETER) $<\
  --mesh=$(word 2, $^)\
  --mesh-hodge-star=$(word 3, $^)\
  --mesh-volumes=$(word 4, $^)\
  --dynamic-library=$(word 2, $|)\
  --input-data=diffusion_transient_advection_diffusion_test_fluid_flow\
  > $@

build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/primal_weak_cochain_brick_3d_2_forman_input.txt:\
  build/$(MODE)/bin/diffusion_transient_discrete_primal_weak_from_continuous$(.EXE)\
  build/$(MODE)/txt/mesh/brick_3d_2/forman.txt\
  build/$(MODE)/txt/mesh/brick_3d_2/forman_vol.txt\
  build/$(MODE)/obj/plugins/diffusion_transient_advection_diffusion_test$(.OBJ)\
  | build/$(MODE)/txt/diffusion/transient/advection_diffusion_test\
    build/$(MODE)/lib/plugins/libdiffusion$(.DLL)
	$(INTERPRETER) $<\
  --mesh=$(word 2, $^)\
  --mesh-volumes=$(word 3, $^)\
  --dynamic-library=$(word 2, $|)\
  --input-data=diffusion_transient_advection_diffusion_test\
  > $@

build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/primal_weak_cochain_brick_3d_2_forman_trapezoidal_0p01_200_potential.txt:\
  build/$(MODE)/bin/advective_diffusion_solve$(.EXE)\
  build/$(MODE)/txt/mesh/brick_3d_2/forman.txt\
  build/$(MODE)/txt/mesh/brick_3d_2/forman_inner.txt\
  build/$(MODE)/txt/mesh/brick_3d_2/forman_inner_corrected.txt\
  build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/mixed_weak_cochain_brick_3d_2_forman_input.txt\
  build/$(MODE)/txt/diffusion/transient/advection_diffusion_test/primal_weak_cochain_brick_3d_2_forman_input.txt\
  | build/$(MODE)/txt/diffusion/transient/advection_diffusion_test
	$(INTERPRETER) $<\
  --mesh=$(word 2, $^)\
  --mesh-inner=$(word 3, $^)\
  --mesh-inner-corrected=$(word 4, $^)\
  --input-data=$(word 5, $^)\
  --input-data-2=$(word 6, $^)\
  --time-step=0.01\
  --number-of-steps=200\
  > $@

txt_advection_diffusion_test:\
  $(_txt_diffusion_transient_advection_diffusion_test)

txt_advection_diffusion_test_clean:
	-$(RM) $(_txt_diffusion_transient_advection_diffusion_test)

txt_advection_diffusion_test_distclean:
	-$(RM) -r build/$(MODE)/txt/diffusion/transient/advection_diffusion_test
  