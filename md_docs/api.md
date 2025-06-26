# api

api

# Application Programming Interface

The PyMOL Application Programming Interface (API) is intended to enable creation of custom software which makes use of PyMOL's molecular visualization, manipulation, and rendering capabilities. Such applications might range from simple Python programs for processing PDB files to advanced molecular modeling suites that integrate multiple computational components with PyMOL. 

An [alphabetical list of cmd methods](/dokuwiki/doku.php?id=api:cmd:alpha "api:cmd:alpha") is available. 

## API Methods vs. Commands

Most PyMOL API methods parallel their PyMOL Command equivalents. For example, the [color](/dokuwiki/doku.php?id=command:color "command:color") command is very similiar in its behavior to the [cmd.color](/dokuwiki/doku.php?id=api:cmd:color "api:cmd:color") API call. The chief differences are: 

\- API methods are called directly from Python whereas commands are parsed via PyMOL's command interpreter. This means that Python programs calling API methods will be dispatched much more quickly than will script-based commands. 

\- The “quiet” argument (if any) is not automatically set to 0 as it is with PyMOL commands. This means that API method output is usually suppressed as compared with the script commands, further increasing efficiency. 

\- Python exceptions may be raised by API methods when things go wrong. 

api.txt · Last modified: 2013/08/19 21:00 (external edit)
  *[API]: Application Programming Interface
