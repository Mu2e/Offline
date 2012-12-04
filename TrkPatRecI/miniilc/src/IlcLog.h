#ifndef ILCLOG_H
#define ILCLOG_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcLog.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

///
/// class for logging debug, info and error messages
///

#include <TObject.h>
#include <TObjArray.h>
#include <TString.h>


class IlcLog: public TObject {
 public:
  IlcLog();
  virtual ~IlcLog();
  static IlcLog* Instance() {return fgInstance;}

  enum EType_t {kFatal = 0, kError, kWarning, kInfo, kDebug, kMaxType};

  static void  EnableDebug(Bool_t enabled);
  static void  SetGlobalLogLevel(EType_t type);
  static Int_t GetGlobalLogLevel();
  static void  SetGlobalDebugLevel(Int_t level);
  static Int_t GetGlobalDebugLevel();
  static void  SetModuleDebugLevel(const char* module, Int_t level);
  static void  ClearModuleDebugLevel(const char* module);
  static void  SetClassDebugLevel(const char* className, Int_t level);
  static void  ClearClassDebugLevel(const char* className);

  static void  SetStandardOutput();
  static void  SetStandardOutput(EType_t type);
  static void  SetErrorOutput();
  static void  SetErrorOutput(EType_t type);
  static void  SetFileOutput(const char* fileName);
  static void  SetFileOutput(EType_t type, const char* fileName);
  static void  Flush();

  static void  SetHandleRootMessages(Bool_t on);

  static void  SetPrintType(Bool_t on);
  static void  SetPrintType(EType_t type, Bool_t on);
  static void  SetPrintModule(Bool_t on);
  static void  SetPrintModule(EType_t type, Bool_t on);
  static void  SetPrintScope(Bool_t on);
  static void  SetPrintScope(EType_t type, Bool_t on);
  static void  SetPrintLocation(Bool_t on);
  static void  SetPrintLocation(EType_t type, Bool_t on);

  static void  SetPrintRepetitions(Bool_t on);

  static void  WriteToFile(const char* name, Int_t option = 0);

  // the following public methods are used by the preprocessor macros 
  // and should not be called directly
  static Bool_t IsDebugEnabled() {return fgDebugEnabled;}
  static Int_t GetDebugLevel(const char* module, const char* className);
  static void  Message(UInt_t level, const char* message, 
                       const char* module, const char* className,
                       const char* function, const char* file, Int_t line);
  static void  Debug(UInt_t level, const char* message, 
                     const char* module, const char* className,
                     const char* function, const char* file, Int_t line);

  static Int_t RedirectStdoutTo(EType_t type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static Int_t RedirectStderrTo(EType_t type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static void  RestoreStdout(Int_t original);
  static void  RestoreStderr(Int_t original);

  static ostream& Stream(EType_t type, UInt_t level,
                         const char* module, const char* className,
                         const char* function, const char* file, Int_t line);

 private:
  IlcLog(const IlcLog& log);
  IlcLog& operator = (const IlcLog& log);

  void           ReadEnvSettings();

  static void    RootErrorHandler(Int_t level, Bool_t abort, 
				  const char* location, const char* message);

  void           CloseFile(Int_t type);
  FILE*          GetOutputStream(Int_t type);

  UInt_t         GetLogLevel(const char* module, const char* className) const;
  void           PrintMessage(UInt_t type, const char* message, 
                              const char* module, const char* className,
                              const char* function, 
                              const char* file, Int_t line);
  void           PrintRepetitions();

  Int_t          RedirectTo(FILE* stream, EType_t type, UInt_t level,
                            const char* module, const char* className,
                            const char* function,
                            const char* file, Int_t line, Bool_t print);

  ostream&       GetStream(EType_t type, UInt_t level,
                           const char* module, const char* className,
                           const char* function, const char* file, Int_t line);

  enum {kDebugOffset = kDebug-1};

  static IlcLog* fgInstance;                 //! pointer to current instance

  static Bool_t  fgDebugEnabled;             // flag for debug en-/disabling

  UInt_t         fGlobalLogLevel;            // global logging level
  TObjArray      fModuleDebugLevels;         // debug levels for modules
  TObjArray      fClassDebugLevels;          // debug levels for classes

  Int_t          fOutputTypes[kMaxType];     // types of output streams
  TString        fFileNames[kMaxType];       // file names
  FILE*          fOutputFiles[kMaxType];     //! log output files
  ofstream*      fOutputStreams[kMaxType];   //! log output streams

  Bool_t         fPrintType[kMaxType];       // print type on/off
  Bool_t         fPrintModule[kMaxType];     // print module on/off
  Bool_t         fPrintScope[kMaxType];      // print scope/class name on/off
  Bool_t         fPrintLocation[kMaxType];   // print file and line on/off

  Bool_t         fPrintRepetitions;          // print number of repetitions instead of repeated message on/off

  Int_t          fRepetitions;               //! counter of repetitions
  UInt_t         fLastType;                  //! type of last message
  TString        fLastMessage;               //! last message
  TString        fLastModule;                //! module name of last message
  TString        fLastClassName;             //! class name of last message
  TString        fLastFunction;              //! function name of last message
  TString        fLastFile;                  //! file name of last message
  Int_t          fLastLine;                  //! line number of last message

  ClassDef(IlcLog, 1)   // class for logging debug, info and error messages
};


// module name
#ifdef __MODULE__
#define MODULENAME() __MODULE__
#else
#define MODULENAME() "NoModule"
#endif

// function name
#if defined(__GNUC__) || defined(__ICC) || defined(__ECC) || defined(__APPLE__)
#define FUNCTIONNAME() __FUNCTION__
// #elif defined(__HP_aCC) || defined(__alpha) || defined(__DECCXX)
// #define FUNCTIONNAME() __FUNC__
#else
#define FUNCTIONNAME() "???"
#endif

// redirection
#define REDIRECTSTDOUT(type, level, scope, whatever) {Int_t originalStdout = IlcLog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; IlcLog::RestoreStdout(originalStdout);}
#define REDIRECTSTDERR(type, level, scope, whatever) {Int_t originalStderr = IlcLog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; IlcLog::RestoreStderr(originalStderr);}
#define REDIRECTSTDOUTANDSTDERR(type, level, scope, whatever) {Int_t originalStdout = IlcLog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); Int_t originalStderr = IlcLog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; IlcLog::RestoreStderr(originalStderr); IlcLog::RestoreStdout(originalStdout);}


// debug level
#ifdef LOG_NO_DEBUG
#define IlcDebugLevel() -1
#define IlcDebugLevelClass() -1
#define IlcDebugLevelGeneral(scope) -1
#else
#define IlcDebugLevel() ((IlcLog::IsDebugEnabled()) ? IlcLog::GetDebugLevel(MODULENAME(), ClassName()) : -1)
#define IlcDebugLevelClass() ((IlcLog::IsDebugEnabled()) ? IlcLog::GetDebugLevel(MODULENAME(), Class()->GetName()) : -1)
#define IlcDebugLevelGeneral(scope) ((IlcLog::IsDebugEnabled()) ? IlcLog::GetDebugLevel(MODULENAME(), scope) : -1)
#endif

// debug messages
#ifdef LOG_NO_DEBUG
#define IlcDebug(level, message)
#define IlcDebugClass(level, message)
#define IlcDebugGeneral(scope, level, message)
#else
#define IlcDebug(level, message) {if (IlcLog::IsDebugEnabled()) IlcLog::Debug(level, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcDebugClass(level, message) {if (IlcLog::IsDebugEnabled()) IlcLog::Debug(level, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcDebugGeneral(scope, level, message) {if (IlcLog::IsDebugEnabled()) IlcLog::Debug(level, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to debug
#define StdoutToIlcDebug(level, whatever) REDIRECTSTDOUT(IlcLog::kDebug, level, ClassName(), whatever)
#define StderrToIlcDebug(level, whatever) REDIRECTSTDERR(IlcLog::kDebug, level, ClassName(), whatever)
#define ToIlcDebug(level, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kDebug, level, ClassName(), whatever)
#define StdoutToIlcDebugClass(level, whatever) REDIRECTSTDOUT(IlcLog::kDebug, level, Class()->GetName(), whatever)
#define StderrToIlcDebugClass(level, whatever) REDIRECTSTDERR(IlcLog::kDebug, level, Class()->GetName(), whatever)
#define ToIlcDebugClass(level, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kDebug, level, Class()->GetName(), whatever)
#define StdoutToIlcDebugGeneral(scope, level, whatever) REDIRECTSTDOUT(IlcLog::kDebug, level, scope, whatever)
#define StderrToIlcDebugGeneral(scope, level, whatever) REDIRECTSTDERR(IlcLog::kDebug, level, scope, whatever)
#define ToIlcDebugGeneral(scope, level, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kDebug, level, scope, whatever)

// debug stream objects
#define IlcDebugStream(level) IlcLog::Stream(IlcLog::kDebug, level, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcDebugClassStream(level) IlcLog::Stream(IlcLog::kDebug, level, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcDebugGeneralStream(scope, level) IlcLog::Stream(IlcLog::kDebug, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// info messages
#ifdef LOG_NO_INFO
#define IlcInfo(message)
#define IlcInfoClass(message)
#define IlcInfoGeneral(scope, message)
#else
#define IlcInfo(message) {IlcLog::Message(IlcLog::kInfo, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcInfoClass(message) {IlcLog::Message(IlcLog::kInfo, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcInfoGeneral(scope, message) {IlcLog::Message(IlcLog::kInfo, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to info
#define StdoutToIlcInfo(whatever) REDIRECTSTDOUT(IlcLog::kInfo, 0, ClassName(), whatever)
#define StderrToIlcInfo(whatever) REDIRECTSTDERR(IlcLog::kInfo, 0, ClassName(), whatever)
#define ToIlcInfo(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kInfo, 0, ClassName(), whatever)
#define StdoutToIlcInfoClass(whatever) REDIRECTSTDOUT(IlcLog::kInfo, 0, Class()->GetName(), whatever)
#define StderrToIlcInfoClass(whatever) REDIRECTSTDERR(IlcLog::kInfo, 0, Class()->GetName(), whatever)
#define ToIlcInfoClass(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kInfo, 0, Class()->GetName(), whatever)
#define StdoutToIlcInfoGeneral(scope, whatever) REDIRECTSTDOUT(IlcLog::kInfo, 0, scope, whatever)
#define StderrToIlcInfoGeneral(scope, whatever) REDIRECTSTDERR(IlcLog::kInfo, 0, scope, whatever)
#define ToIlcInfoGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kInfo, 0, scope, whatever)

// info stream objects
#define IlcInfoStream() IlcLog::Stream(IlcLog::kInfo, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcInfoClassStream() IlcLog::Stream(IlcLog::kInfo, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcInfoGeneralStream(scope) IlcLog::Stream(IlcLog::kInfo, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// warning messages
#ifdef LOG_NO_WARNING
#define IlcWarning(message)
#define IlcWarningClass(message)
#define IlcWarningGeneral(scope, message)
#else
#define IlcWarning(message) {IlcLog::Message(IlcLog::kWarning, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcWarningClass(message) {IlcLog::Message(IlcLog::kWarning, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcWarningGeneral(scope, message) {IlcLog::Message(IlcLog::kWarning, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to warning
#define StdoutToIlcWarning(whatever) REDIRECTSTDOUT(IlcLog::kWarning, 0, ClassName(), whatever)
#define StderrToIlcWarning(whatever) REDIRECTSTDERR(IlcLog::kWarning, 0, ClassName(), whatever)
#define ToIlcWarning(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kWarning, 0, ClassName(), whatever)
#define StdoutToIlcWarningClass(whatever) REDIRECTSTDOUT(IlcLog::kWarning, 0, Class()->GetName(), whatever)
#define StderrToIlcWarningClass(whatever) REDIRECTSTDERR(IlcLog::kWarning, 0, Class()->GetName(), whatever)
#define ToIlcWarningClass(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kWarning, 0, Class()->GetName(), whatever)
#define StdoutToIlcWarningGeneral(scope, whatever) REDIRECTSTDOUT(IlcLog::kWarning, 0, scope, whatever)
#define StderrToIlcWarningGeneral(scope, whatever) REDIRECTSTDERR(IlcLog::kWarning, 0, scope, whatever)
#define ToIlcWarningGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kWarning, 0, scope, whatever)

// warning stream objects
#define IlcWarningStream() IlcLog::Stream(IlcLog::kWarning, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcWarningClassStream() IlcLog::Stream(IlcLog::kWarning, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcWarningGeneralStream(scope) IlcLog::Stream(IlcLog::kWarning, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// error messages
#define IlcError(message) {IlcLog::Message(IlcLog::kError, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcErrorClass(message) {IlcLog::Message(IlcLog::kError, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcErrorGeneral(scope, message) {IlcLog::Message(IlcLog::kError, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}

// redirection to error
#define StdoutToIlcError(whatever) REDIRECTSTDOUT(IlcLog::kError, 0, ClassName(), whatever)
#define StderrToIlcError(whatever) REDIRECTSTDERR(IlcLog::kError, 0, ClassName(), whatever)
#define ToIlcError(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kError, 0, ClassName(), whatever)
#define StdoutToIlcErrorClass(whatever) REDIRECTSTDOUT(IlcLog::kError, 0, Class()->GetName(), whatever)
#define StderrToIlcErrorClass(whatever) REDIRECTSTDERR(IlcLog::kError, 0, Class()->GetName(), whatever)
#define ToIlcErrorClass(whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kError, 0, Class()->GetName(), whatever)
#define StdoutToIlcErrorGeneral(scope, whatever) REDIRECTSTDOUT(IlcLog::kError, 0, scope, whatever)
#define StderrToIlcErrorGeneral(scope, whatever) REDIRECTSTDERR(IlcLog::kError, 0, scope, whatever)
#define ToIlcErrorGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(IlcLog::kError, 0, scope, whatever)

// error stream objects
#define IlcErrorStream() IlcLog::Stream(IlcLog::kError, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcErrorClassStream() IlcLog::Stream(IlcLog::kError, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define IlcErrorGeneralStream(scope) IlcLog::Stream(IlcLog::kError, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// fatal messages
#define IlcFatal(message) {IlcLog::Message(IlcLog::kFatal, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcFatalClass(message) {IlcLog::Message(IlcLog::kFatal, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define IlcFatalGeneral(scope, message) {IlcLog::Message(IlcLog::kFatal, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}

#endif
