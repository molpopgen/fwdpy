TODO 
=================

* Document custom temporal samplers.  Target version 0.0.4
* Add safety check to (de)serialization.  Current method trusts the input string completely.  This likely requires a
  saftey-checking template that is specialized for each population type.  Target version 0.0.5.
* Remove fwdpy.fwdpyio module in favor of object-level serialization.  Target version 0.0.5.
* Serialization should be at the level of PopType and PopVec objects.  The latter can be done in parallel.  Target
  version 0.0.5.
* Serialization should support direct to/from file, with gzip as option via zlib. Target version 0.0.5.
* Custom rules classes. This will allow "stateful" fitness models, and many other things.  Target version 0.0.5.
