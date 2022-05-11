// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME ElissaTreeEventDict
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
#include "ElissaTreeEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_ElissaTreeEvent(void *p = 0);
   static void *newArray_ElissaTreeEvent(Long_t size, void *p);
   static void delete_ElissaTreeEvent(void *p);
   static void deleteArray_ElissaTreeEvent(void *p);
   static void destruct_ElissaTreeEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ElissaTreeEvent*)
   {
      ::ElissaTreeEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ElissaTreeEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ElissaTreeEvent", ::ElissaTreeEvent::Class_Version(), "ElissaTreeEvent.h", 10,
                  typeid(::ElissaTreeEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ElissaTreeEvent::Dictionary, isa_proxy, 4,
                  sizeof(::ElissaTreeEvent) );
      instance.SetNew(&new_ElissaTreeEvent);
      instance.SetNewArray(&newArray_ElissaTreeEvent);
      instance.SetDelete(&delete_ElissaTreeEvent);
      instance.SetDeleteArray(&deleteArray_ElissaTreeEvent);
      instance.SetDestructor(&destruct_ElissaTreeEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ElissaTreeEvent*)
   {
      return GenerateInitInstanceLocal((::ElissaTreeEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ElissaTreeEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ElissaTreeEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ElissaTreeEvent::Class_Name()
{
   return "ElissaTreeEvent";
}

//______________________________________________________________________________
const char *ElissaTreeEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ElissaTreeEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ElissaTreeEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ElissaTreeEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ElissaTreeEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ElissaTreeEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ElissaTreeEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ElissaTreeEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ElissaTreeEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class ElissaTreeEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ElissaTreeEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(ElissaTreeEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ElissaTreeEvent(void *p) {
      return  p ? new(p) ::ElissaTreeEvent : new ::ElissaTreeEvent;
   }
   static void *newArray_ElissaTreeEvent(Long_t nElements, void *p) {
      return p ? new(p) ::ElissaTreeEvent[nElements] : new ::ElissaTreeEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_ElissaTreeEvent(void *p) {
      delete ((::ElissaTreeEvent*)p);
   }
   static void deleteArray_ElissaTreeEvent(void *p) {
      delete [] ((::ElissaTreeEvent*)p);
   }
   static void destruct_ElissaTreeEvent(void *p) {
      typedef ::ElissaTreeEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ElissaTreeEvent

namespace ROOT {
   static TClass *vectorlEElissaTreeEventgR_Dictionary();
   static void vectorlEElissaTreeEventgR_TClassManip(TClass*);
   static void *new_vectorlEElissaTreeEventgR(void *p = 0);
   static void *newArray_vectorlEElissaTreeEventgR(Long_t size, void *p);
   static void delete_vectorlEElissaTreeEventgR(void *p);
   static void deleteArray_vectorlEElissaTreeEventgR(void *p);
   static void destruct_vectorlEElissaTreeEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ElissaTreeEvent>*)
   {
      vector<ElissaTreeEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ElissaTreeEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ElissaTreeEvent>", -2, "vector", 386,
                  typeid(vector<ElissaTreeEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEElissaTreeEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ElissaTreeEvent>) );
      instance.SetNew(&new_vectorlEElissaTreeEventgR);
      instance.SetNewArray(&newArray_vectorlEElissaTreeEventgR);
      instance.SetDelete(&delete_vectorlEElissaTreeEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEElissaTreeEventgR);
      instance.SetDestructor(&destruct_vectorlEElissaTreeEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ElissaTreeEvent> >()));

      ::ROOT::AddClassAlternate("vector<ElissaTreeEvent>","std::vector<ElissaTreeEvent, std::allocator<ElissaTreeEvent> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ElissaTreeEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEElissaTreeEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ElissaTreeEvent>*)0x0)->GetClass();
      vectorlEElissaTreeEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEElissaTreeEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEElissaTreeEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ElissaTreeEvent> : new vector<ElissaTreeEvent>;
   }
   static void *newArray_vectorlEElissaTreeEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ElissaTreeEvent>[nElements] : new vector<ElissaTreeEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEElissaTreeEventgR(void *p) {
      delete ((vector<ElissaTreeEvent>*)p);
   }
   static void deleteArray_vectorlEElissaTreeEventgR(void *p) {
      delete [] ((vector<ElissaTreeEvent>*)p);
   }
   static void destruct_vectorlEElissaTreeEventgR(void *p) {
      typedef vector<ElissaTreeEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ElissaTreeEvent>

namespace {
  void TriggerDictionaryInitialization_ElissaTreeEventDict_Impl() {
    static const char* headers[] = {
"ElissaTreeEvent.h",
0
    };
    static const char* includePaths[] = {
"/home/testov/root/include/",
"/home/testov/DelilaSorting/lib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "ElissaTreeEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ElissaTreeEvent.h")))  ElissaTreeEvent;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "ElissaTreeEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "ElissaTreeEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ElissaTreeEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ElissaTreeEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ElissaTreeEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ElissaTreeEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ElissaTreeEventDict() {
  TriggerDictionaryInitialization_ElissaTreeEventDict_Impl();
}
