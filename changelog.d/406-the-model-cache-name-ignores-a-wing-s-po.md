### Fixed

- The model bin name and its validity hash now carry the concrete types of the wings' aero models and the winches' models, so two builds whose polar interpolants (e.g. `remove_nan`, polar format, panel count) or winch model differ each get their own bin instead of the second failing in `write_callable!` on the first's.
