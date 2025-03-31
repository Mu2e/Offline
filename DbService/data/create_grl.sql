-- create tables for goodrun list in the dqm database
-- should run as admin_role, which should own all objects
SET ROLE admin_role;

--
--
CREATE SCHEMA grl;
GRANT USAGE ON SCHEMA grl TO PUBLIC;

-- definitons of a grl word
CREATE TABLE grl.words 
  (name TEXT NOT NULL,
  description TEXT,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT words_pk PRIMARY KEY (name) );
GRANT SELECT ON grl.words TO PUBLIC;
GRANT INSERT ON grl.words TO grlwrite;

-- definitons of a grl bits within a word
CREATE TABLE grl.bits 
  (name TEXT NOT NULL,
  bitname TEXT NOT NULL,
  bitnumber INTEGER NOT NULL,
  description TEXT,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT bits_pk PRIMARY KEY (name,bitname),
  CONSTRAINT bits_name_fk FOREIGN KEY (name) REFERENCES grl.words(name));
  GRANT SELECT ON grl.bits TO PUBLIC;
  GRANT INSERT ON grl.bits TO grlwrite;

-- grl entries:  word, value and IoV
CREATE TABLE grl.entries 
  (eid SERIAL,
  name TEXT NOT NULL,
  iov TEXT NOT NULL,
  value INTEGER NOT NULL,
  retired INTEGER NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT entries_name_fk FOREIGN KEY (name) REFERENCES grl.words(name));
  GRANT SELECT ON grl.entries TO PUBLIC;
  GRANT INSERT ON grl.entries TO grlwrite;
  GRANT UPDATE ON grl.entries TO grlwrite;
  GRANT UPDATE ON grl.entries_eid_seq TO grlwrite;

-- grl list definitions
CREATE TABLE grl.lists
  (lid SERIAL,
  name TEXT NOT NULL,
  locked INTEGER NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT lists_pk PRIMARY KEY (name),
  CONSTRAINT lists_lid_uc UNIQUE (lid) );
  GRANT SELECT ON grl.lists TO PUBLIC;
  GRANT INSERT ON grl.lists TO grlwrite;
  GRANT UPDATE ON grl.lists TO grlwrite;
  GRANT UPDATE ON grl.lists_lid_seq TO grlwrite;

-- grl list entries
CREATE TABLE grl.iovs
  (lid INTEGER NOT NULL,
  iov TEXT NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,
  CONSTRAINT iovs_lid_fk FOREIGN KEY (lid) REFERENCES grl.lists(lid) );
  GRANT SELECT ON grl.iovs TO PUBLIC;
  GRANT INSERT ON grl.iovs TO grlwrite;

