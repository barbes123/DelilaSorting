// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LaBrTreeEventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LaBrTreeEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LaBrTreeEvent(void *p = 0);
   static void *newArray_LaBrTreeEvent(Long_t size, void *p);
   static void delete_LaBrTreeEvent(void *p);
   static void deleteArray_LaBrTreeEvent(void *p);
   static void destruct_LaBrTreeEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LaBrTreeEvent*)
   {
      ::LaBrTreeEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LaBrTreeEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LaBrTreeEvent", ::LaBrTreeEvent::Class_Version(), "LaBrTreeEvent.h", 10,
                  typeid(::LaBrTreeEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LaBrTreeEvent::Dictionary, isa_proxy, 4,
                  sizeof(::LaBrTreeEvent) );
      instance.SetNew(&new_LaBrTreeEvent);
      instance.SetNewArray(&newArray_LaBrTreeEvent);
      instance.SetDelete(&delete_LaBrTreeEvent);
      instance.SetDeleteArray(&deleteArray_LaBrTreeEvent);
      instance.SetDestructor(&destruct_LaBrTreeEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LaBrTreeEvent*)
   {
      return GenerateInitInstanceLocal((::LaBrTreeEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LaBrTreeEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LaBrTreeEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LaBrTreeEvent::Class_Name()
{
   return "LaBrTreeEvent";
}

//______________________________________________________________________________
const char *LaBrTreeEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LaBrTreeEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LaBrTreeEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LaBrTreeEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LaBrTreeEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LaBrTreeEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LaBrTreeEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LaBrTreeEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LaBrTreeEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class LaBrTreeEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LaBrTreeEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(LaBrTreeEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LaBrTreeEvent(void *p) {
      return  p ? new(p) ::LaBrTreeEvent : new ::LaBrTreeEvent;
   }
   static void *newArray_LaBrTreeEvent(Long_t nElements, void *p) {
      return p ? new(p) ::LaBrTreeEvent[nElements] : new ::LaBrTreeEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_LaBrTreeEvent(void *p) {
      delete ((::LaBrTreeEvent*)p);
   }
   static void deleteArray_LaBrTreeEvent(void *p) {
      delete [] ((::LaBrTreeEvent*)p);
   }
   static void destruct_LaBrTreeEvent(void *p) {
      typedef ::LaBrTreeEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LaBrTreeEvent

namespace ROOT {
   static TClass *vectorlELaBrTreeEventgR_Dictionary();
   static void vectorlELaBrTreeEventgR_TClassManip(TClass*);
   static void *new_vectorlELaBrTreeEventgR(void *p = 0);
   static void *newArray_vectorlELaBrTreeEventgR(Long_t size, void *p);
   static void delete_vectorlELaBrTreeEventgR(void *p);
   static void deleteArray_vectorlELaBrTreeEventgR(void *p);
   static void destruct_vectorlELaBrTreeEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LaBrTreeEvent>*)
   {
      vector<LaBrTreeEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LaBrTreeEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LaBrTreeEvent>", -2, "vector", 386,
                  typeid(vector<LaBrTreeEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELaBrTreeEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LaBrTreeEvent>) );
      instance.SetNew(&new_vectorlELaBrTreeEventgR);
      instance.SetNewArray(&newArray_vectorlELaBrTreeEventgR);
      instance.SetDelete(&delete_vectorlELaBrTreeEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlELaBrTreeEventgR);
      instance.SetDestructor(&destruct_vectorlELaBrTreeEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LaBrTreeEvent> >()));

      ::ROOT::AddClassAlternate("vector<LaBrTreeEvent>","std::vector<LaBrTreeEvent, std::allocator<LaBrTreeEvent> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LaBrTreeEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELaBrTreeEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LaBrTreeEvent>*)0x0)->GetClass();
      vectorlELaBrTreeEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELaBrTreeEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELaBrTreeEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LaBrTreeEvent> : new vector<LaBrTreeEvent>;
   }
   static void *newArray_vectorlELaBrTreeEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LaBrTreeEvent>[nElements] : new vector<LaBrTreeEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELaBrTreeEventgR(void *p) {
      delete ((vector<LaBrTreeEvent>*)p);
   }
   static void deleteArray_vectorlELaBrTreeEventgR(void *p) {
      delete [] ((vector<LaBrTreeEvent>*)p);
   }
   static void destruct_vectorlELaBrTreeEventgR(void *p) {
      typedef vector<LaBrTreeEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LaBrTreeEvent>

namespace {
  void TriggerDictionaryInitialization_LaBrTreeEventDict_Impl() {
    static const char* headers[] = {
"LaBrTreeEvent.h",
0
    };
    static const char* includePaths[] = {
"/home/testov/root/include/",
"/home/testov/DelilaSorting/lib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LaBrTreeEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$LaBrTreeEvent.h")))  LaBrTreeEvent;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LaBrTreeEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "LaBrTreeEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"LaBrTreeEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LaBrTreeEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LaBrTreeEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LaBrTreeEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LaBrTreeEventDict() {
  TriggerDictionaryInitialization_LaBrTreeEventDict_Impl();
}
