// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DelilaEventDict
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
#include "DelilaEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *DelilaEvent_Dictionary();
   static void DelilaEvent_TClassManip(TClass*);
   static void *new_DelilaEvent(void *p = 0);
   static void *newArray_DelilaEvent(Long_t size, void *p);
   static void delete_DelilaEvent(void *p);
   static void deleteArray_DelilaEvent(void *p);
   static void destruct_DelilaEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DelilaEvent*)
   {
      ::DelilaEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DelilaEvent));
      static ::ROOT::TGenericClassInfo 
         instance("DelilaEvent", "DelilaEvent.h", 10,
                  typeid(::DelilaEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DelilaEvent_Dictionary, isa_proxy, 4,
                  sizeof(::DelilaEvent) );
      instance.SetNew(&new_DelilaEvent);
      instance.SetNewArray(&newArray_DelilaEvent);
      instance.SetDelete(&delete_DelilaEvent);
      instance.SetDeleteArray(&deleteArray_DelilaEvent);
      instance.SetDestructor(&destruct_DelilaEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DelilaEvent*)
   {
      return GenerateInitInstanceLocal((::DelilaEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DelilaEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DelilaEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DelilaEvent*)0x0)->GetClass();
      DelilaEvent_TClassManip(theClass);
   return theClass;
   }

   static void DelilaEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_DelilaEvent(void *p) {
      return  p ? new(p) ::DelilaEvent : new ::DelilaEvent;
   }
   static void *newArray_DelilaEvent(Long_t nElements, void *p) {
      return p ? new(p) ::DelilaEvent[nElements] : new ::DelilaEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_DelilaEvent(void *p) {
      delete ((::DelilaEvent*)p);
   }
   static void deleteArray_DelilaEvent(void *p) {
      delete [] ((::DelilaEvent*)p);
   }
   static void destruct_DelilaEvent(void *p) {
      typedef ::DelilaEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DelilaEvent

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 386,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 386,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEDelilaEventgR_Dictionary();
   static void vectorlEDelilaEventgR_TClassManip(TClass*);
   static void *new_vectorlEDelilaEventgR(void *p = 0);
   static void *newArray_vectorlEDelilaEventgR(Long_t size, void *p);
   static void delete_vectorlEDelilaEventgR(void *p);
   static void deleteArray_vectorlEDelilaEventgR(void *p);
   static void destruct_vectorlEDelilaEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DelilaEvent>*)
   {
      vector<DelilaEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DelilaEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DelilaEvent>", -2, "vector", 386,
                  typeid(vector<DelilaEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDelilaEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<DelilaEvent>) );
      instance.SetNew(&new_vectorlEDelilaEventgR);
      instance.SetNewArray(&newArray_vectorlEDelilaEventgR);
      instance.SetDelete(&delete_vectorlEDelilaEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDelilaEventgR);
      instance.SetDestructor(&destruct_vectorlEDelilaEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DelilaEvent> >()));

      ::ROOT::AddClassAlternate("vector<DelilaEvent>","std::vector<DelilaEvent, std::allocator<DelilaEvent> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DelilaEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDelilaEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DelilaEvent>*)0x0)->GetClass();
      vectorlEDelilaEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDelilaEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDelilaEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DelilaEvent> : new vector<DelilaEvent>;
   }
   static void *newArray_vectorlEDelilaEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DelilaEvent>[nElements] : new vector<DelilaEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDelilaEventgR(void *p) {
      delete ((vector<DelilaEvent>*)p);
   }
   static void deleteArray_vectorlEDelilaEventgR(void *p) {
      delete [] ((vector<DelilaEvent>*)p);
   }
   static void destruct_vectorlEDelilaEventgR(void *p) {
      typedef vector<DelilaEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DelilaEvent>

namespace {
  void TriggerDictionaryInitialization_DelilaEventDict_Impl() {
    static const char* headers[] = {
"DelilaEvent.h",
0
    };
    static const char* includePaths[] = {
"/home/testov/root/include/",
"/home/testov/DelilaSorting/lib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DelilaEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$DelilaEvent.h")))  DelilaEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DelilaEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "DelilaEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"DelilaEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DelilaEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DelilaEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DelilaEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DelilaEventDict() {
  TriggerDictionaryInitialization_DelilaEventDict_Impl();
}
