### Added
- `save_log(logger, sys, name)` saves a log that carries the structure document of
  `sys` as JSON under the metadata key `topology`, uncompressed so that arrow-js can
  read it in a browser. `load_log` returns the document in `SysLog.metadata`.
